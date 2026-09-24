% Runs every case in a cases JSON file through CUFSM's own MATLAB code and writes what it
% computes, stage by stage, as JSON for the Rust tests to compare against.
%
%   octave-cli --no-gui run_cases.m <cufsm-root> <cases.json> <out.json>
%
% <cufsm-root> is a checkout of CUFSM (https://github.com/thinwalled/cufsm-git). Nothing in it is
% modified: the shims/ folder stands in for the progress bar stripmain.m draws.

args = argv();
root = args{1};
cases_file = args{2};
out_file = args{3};
here = fileparts(mfilename('fullpath'));
% Octave's own eigs.m, renamed, for the eigs shim to call. It is GPL, so it is copied from the
% local Octave install at run time into a temporary folder and never kept in this repository.
octave_eigs = which('eigs');
tmp = tempname(); mkdir(tmp);
src = fileread(octave_eigs);
src = regexprep(src, 'function varargout = eigs \(', 'function varargout = eigs_octave (', 'once');
fid = fopen(fullfile(tmp, 'eigs_octave.m'), 'w'); fputs(fid, src); fclose(fid);
addpath(tmp);
addpath(fullfile(here, 'shims'));
addpath(fullfile(root, 'analysis'));
addpath(fullfile(root, 'analysis', 'cFSM'));
addpath(fullfile(root, 'helpers'));

cases = jsondecode(fileread(cases_file));
if ~iscell(cases)
  cases = num2cell(cases);
end
out = cell(numel(cases), 1);
for c = 1:numel(cases)
  cs = cases{c};
  E = cs.E; nu = cs.nu;
  if isfield(cs, 'template')
    tp = cs.template;
    [prop, node, elem] = templatecalc(tp.CorZ, tp.h, tp.b1, tp.b2, tp.d1, tp.d2, tp.r1, tp.r2, tp.r3, tp.r4, ...
      tp.q1, tp.q2, tp.t, tp.nh, tp.nb1, tp.nb2, tp.nd1, tp.nd2, tp.nr1, tp.nr2, tp.nr3, tp.nr4, tp.kipin, tp.center);
    prop = [100 E E nu nu E / (2 * (1 + nu))];
    elem(:, 5) = 100;
  else
    node = cs.node;
    elem = cs.elem;
    prop = [100 E E nu nu E / (2 * (1 + nu))];
    elem(:, 5) = 100;
  end
  % Fixed DOFs: rows [node# column], column 4..7 = dofx dofz dofy dofrot.
  if isfield(cs, 'fix') && ~isempty(cs.fix)
    fx = cs.fix;
    if isvector(fx) && numel(fx) == 2, fx = fx(:)'; end
    for i = 1:size(fx, 1)
      node(fx(i, 1), fx(i, 2)) = 0;
    end
  end
  springs = 0;
  if isfield(cs, 'springs') && ~isempty(cs.springs)
    springs = cs.springs;
    if isvector(springs), springs = springs(:)'; end
  end
  constraints = 0;
  if isfield(cs, 'constraints') && ~isempty(cs.constraints)
    constraints = cs.constraints;
    if isvector(constraints), constraints = constraints(:)'; end
  end
  % Reference stresses from the actions, as the CUFSM GUI does.
  [A, xcg, zcg, Ixx, Izz, Ixz, thetap, I11, I22] = grosprop(node, elem);
  ac = cs.actions;
  node = stresgen(node, ac.P, ac.Mxx, ac.Mzz, ac.M11, ac.M22, A, xcg, zcg, Ixx, Izz, Ixz, thetap, I11, I22, ac.unsymm);

  lengths = cs.lengths(:)';
  if iscell(cs.m_all)
    m_all = cellfun(@(v) v(:)', cs.m_all, 'UniformOutput', false)';
  else
    mm = cs.m_all;
    if isvector(mm) && numel(lengths) > 1 && numel(mm) == numel(lengths), mm = mm(:); end
    m_all = cell(1, numel(lengths));
    for l = 1:numel(lengths), m_all{l} = mm(l, :); end
  end
  BC = cs.bc;
  neigs = cs.neigs;
  GBTcon.glob = 0; GBTcon.dist = 0; GBTcon.local = 0; GBTcon.other = 0;
  GBTcon.ospace = 1; GBTcon.couple = 1; GBTcon.orth = 1; GBTcon.norm = 0;

  % Stage 1: strip 1's local matrices and their global rotation, at the first length.
  a = lengths(1);
  m_a = msort({m_all{1}}){1};
  elprop = elemprop(node, elem, size(node, 1), size(elem, 1));
  t = elem(1, 4); b = elprop(1, 2);
  kl = klocal(E, E, nu, nu, E / (2 * (1 + nu)), t, a, b, BC, m_a);
  kgl = kglocal(a, b, node(elem(1, 2), 8) * t, node(elem(1, 3), 8) * t, BC, m_a);
  [kgk, kgkg] = trans(elprop(1, 3), kl, kgl, m_a);

  % Stage 2: the assembled global matrices at the first length (small models only).
  nn = size(node, 1); tm = numel(m_a);
  K = []; Kg = [];
  if 4 * nn * tm <= 400
    K = sparse(zeros(4 * nn * tm)); Kg = sparse(zeros(4 * nn * tm));
    for i = 1:size(elem, 1)
      ti = elem(i, 4); bi = elprop(i, 2);
      k_l = klocal(E, E, nu, nu, E / (2 * (1 + nu)), ti, a, bi, BC, m_a);
      kg_l = kglocal(a, bi, node(elem(i, 2), 8) * ti, node(elem(i, 3), 8) * ti, BC, m_a);
      [k, kg] = trans(elprop(i, 3), k_l, kg_l, m_a);
      [K, Kg] = assemble(K, Kg, k, kg, elem(i, 2), elem(i, 3), nn, m_a);
    end
    K = full(K); Kg = full(Kg);
  end

  % Stage 3: the full analysis, CUFSM's stripmain.m unmodified.
  [curve, shapes] = stripmain(prop, node, elem, lengths, springs, constraints, GBTcon, BC, m_all, neigs);
  lf = cell(numel(lengths), 1);
  mode1 = cell(numel(lengths), 1);
  for l = 1:numel(lengths)
    v = curve{l}(:, 2);
    lf{l} = v(v > 0)';
    s = shapes{l};
    if isempty(s), mode1{l} = []; else, mode1{l} = s(:, 1)'; end
  end

  % Stage 3b: the same reduced problem solved with the full eig() instead of eigs(): CUFSM's own
  % assembly (as stripmain.m) and constraint basis (constr_user + null, as stripmain.m), so the
  % load factors come without ARPACK's convergence tolerance. The commented-out solver path in
  % stripmain.m used eig() the same way.
  lf_dense = cell(numel(lengths), 1);
  for l = 1:numel(lengths)
    al = lengths(l);
    ma = msort({m_all{l}}){1};
    tml = numel(ma);
    KL = sparse(zeros(4 * nn * tml)); KgL = sparse(zeros(4 * nn * tml));
    for i = 1:size(elem, 1)
      ti = elem(i, 4); bi = elprop(i, 2);
      k_l = klocal(E, E, nu, nu, E / (2 * (1 + nu)), ti, al, bi, BC, ma);
      kg_l = kglocal(al, bi, node(elem(i, 2), 8) * ti, node(elem(i, 3), 8) * ti, BC, ma);
      [k, kg] = trans(elprop(i, 3), k_l, kg_l, ma);
      [KL, KgL] = assemble(KL, KgL, k, kg, elem(i, 2), elem(i, 3), nn, ma);
    end
    % Springs exactly as stripmain.m adds them (the v4.3 method).
    if ~isempty(springs) && size(springs, 2) == 10 && springs(1, 1) ~= 0
      for si = 1:size(springs, 1)
        ks_l = spring_klocal(springs(si, 4), springs(si, 5), springs(si, 6), springs(si, 7), al, BC, ma, springs(si, 9), springs(si, 10) * al);
        ni_s = springs(si, 2); nj_s = springs(si, 3);
        if nj_s == 0
          alpha_s = 0;
        else
          dxs = node(nj_s, 2) - node(ni_s, 2); dzs = node(nj_s, 3) - node(ni_s, 3);
          if sqrt(dxs^2 + dzs^2) < 1e-10 || springs(si, 8) == 0
            alpha_s = 0;
          else
            alpha_s = atan2(dzs, dxs);
          end
        end
        ks = spring_trans(alpha_s, ks_l, ma);
        KL = spring_assemble(KL, ks, ni_s, nj_s, nn, ma);
      end
    end
    if constr_BCFlag(node, constraints) == 0
      R = speye(4 * nn * tml);
    else
      Ruser = constr_user(node, constraints, ma);
      Ru0 = null(Ruser');
      R = null(Ru0');
    end
    Kff = full(R' * KL * R); Kgff = full(R' * KgL * R);
    ev = eig(Kff, (Kgff + Kgff') / 2);
    ev = real(ev(abs(imag(ev)) < 1e-5 * abs(real(ev)) & real(ev) > 0));
    ev = sort(ev);
    lf_dense{l} = ev(1:min(neigs, numel(ev)))';
  end

  % Stage 4 (cases flagged cfsm): CUFSM's cFSM. The natural basis at the first length and its
  % space sizes; load factors restricted to G, D and L alone at every length (CUFSM's base vectors
  % and mode_select, solved with eig()); and the default classification (ospace 1, couple 1,
  % orth 2, norm 1) of the unconstrained modes of the dense solve.
  cf = struct();
  if isfield(cs, 'cfsm') && cs.cfsm
    [cutA, cutxc, cutzc, cutIx, cutIz, cutIxz, cuttheta, cutI1, cutI2, cutJ, cutxs, cutzs, cutCw, cutB1, cutB2, cutwn] = cutwp_prop2(node(:, 2:3), elem(:, 2:4));
    cf.cutwp = struct('A', cutA, 'xc', cutxc, 'zc', cutzc, 'Ix', cutIx, 'Iz', cutIz, 'Ixz', cutIxz, 'theta', cuttheta, 'I1', cutI1, 'I2', cutI2, 'J', cutJ, 'xs', cutxs, 'zs', cutzs, 'Cw', cutCw, 'B1', cutB1, 'B2', cutB2, 'wn', cutwn(:)');
    ma1 = msort({m_all{1}}){1};
    [bvl, ngm, ndm, nlm] = base_column(node, elem, prop, lengths(1), BC, ma1);
    cf.ngm = ngm; cf.ndm = ndm; cf.nlm = nlm;
    cf.b_v_l = full(bvl);
    spaces = {'G', [1 0 0 0]; 'D', [0 1 0 0]; 'L', [0 0 1 0]};
    for sp = 1:size(spaces, 1)
      lfs = cell(numel(lengths), 1);
      for l = 1:numel(lengths)
        al = lengths(l); ma = msort({m_all{l}}){1};
        [bv, ng, nd, nl] = base_column(node, elem, prop, al, BC, ma);
        fl = spaces{sp, 2};
        % CUFSM's mode_select fails on an empty space (b_v_red_m never assigned): a plain channel
        % has no distortional modes. Recorded as no load factors.
        if [ng nd nl] * fl(1:3)' == 0, lfs{l} = []; continue; end
        R = mode_select(bv, ng, nd, nl, ones(1, ng) * fl(1), ones(1, nd) * fl(2), ones(1, nl) * fl(3), ones(1, 4 * nn - ng - nd - nl) * fl(4), 4 * nn, ma);
        KL = sparse(zeros(4 * nn * numel(ma))); KgL = KL;
        for i = 1:size(elem, 1)
          ti = elem(i, 4); bi = elprop(i, 2);
          k_l = klocal(E, E, nu, nu, E / (2 * (1 + nu)), ti, al, bi, BC, ma);
          kg_l = kglocal(al, bi, node(elem(i, 2), 8) * ti, node(elem(i, 3), 8) * ti, BC, ma);
          [k, kg] = trans(elprop(i, 3), k_l, kg_l, ma);
          [KL, KgL] = assemble(KL, KgL, k, kg, elem(i, 2), elem(i, 3), nn, ma);
        end
        % Springs into K and constraints into R, as stripmain.m does with modal constraints on.
        if ~isempty(springs) && size(springs, 2) == 10 && springs(1, 1) ~= 0
          for si = 1:size(springs, 1)
            ks_l = spring_klocal(springs(si, 4), springs(si, 5), springs(si, 6), springs(si, 7), al, BC, ma, springs(si, 9), springs(si, 10) * al);
            ni_s = springs(si, 2); nj_s = springs(si, 3);
            alpha_s = 0;
            if nj_s ~= 0
              dxs = node(nj_s, 2) - node(ni_s, 2); dzs = node(nj_s, 3) - node(ni_s, 3);
              if sqrt(dxs^2 + dzs^2) >= 1e-10 && springs(si, 8) ~= 0, alpha_s = atan2(dzs, dxs); end
            end
            KL = spring_assemble(KL, spring_trans(alpha_s, ks_l, ma), ni_s, nj_s, nn, ma);
          end
        end
        if constr_BCFlag(node, constraints) ~= 0
          Rm0 = null(R');
          Ru0 = null(constr_user(node, constraints, ma)');
          R = null([Rm0 Ru0]');
        end
        if isempty(R) || size(R, 2) == 0, lfs{l} = []; continue; end
        ev = eig(full(R' * KL * R), full((R' * KgL * R + (R' * KgL * R)') / 2));
        ev = real(ev(abs(imag(ev)) < 1e-5 * abs(real(ev)) & real(ev) > 0));
        ev = sort(ev);
        lfs{l} = ev(1:min(neigs, numel(ev)))';
      end
      cf.(['lf_' spaces{sp, 1}]) = lfs;
    end
    % Classification of the unconstrained modes (dense solve), CUFSM's defaults.
    G2.ospace = 1; G2.couple = 1; G2.orth = 2; G2.norm = 1;
    clas = cell(numel(lengths), 1);
    clas_nat = cell(numel(lengths), 1);
    for l = 1:numel(lengths)
      al = lengths(l); ma = msort({m_all{l}}){1};
      KL = sparse(zeros(4 * nn * numel(ma))); KgL = KL;
      for i = 1:size(elem, 1)
        ti = elem(i, 4); bi = elprop(i, 2);
        k_l = klocal(E, E, nu, nu, E / (2 * (1 + nu)), ti, al, bi, BC, ma);
        kg_l = kglocal(al, bi, node(elem(i, 2), 8) * ti, node(elem(i, 3), 8) * ti, BC, ma);
        [k, kg] = trans(elprop(i, 3), k_l, kg_l, ma);
        [KL, KgL] = assemble(KL, KgL, k, kg, elem(i, 2), elem(i, 3), nn, ma);
      end
      [V, D] = eig(full(KL), full((KgL + KgL') / 2));
      ev = diag(D); ok = find(abs(imag(ev)) < 1e-5 * abs(real(ev)) & real(ev) > 0);
      [~, ord] = sort(real(ev(ok))); ok = ok(ord); ok = ok(1:min(3, numel(ok)));
      [bv, ng, nd, nl] = base_column(node, elem, prop, al, BC, ma);
      bvu = base_update(G2.ospace, G2.norm, bv, al, ma, node, elem, prop, ng, nd, nl, BC, G2.couple, G2.orth);
      clq = zeros(numel(ok), 4);
      for q = 1:numel(ok)
        clq(q, :) = mode_class(bvu, real(V(:, ok(q))), ng, nd, nl, ma, 4 * nn, G2.couple);
      end
      clas{l} = clq;
      % The same with the natural basis (orth 1): no eig() in base_update, so the basis is unique
      % wherever the distortional space is (ndm <= 1, or none).
      bvn = base_update(1, 1, bv, al, ma, node, elem, prop, ng, nd, nl, BC, 1, 1);
      cln = zeros(numel(ok), 4);
      for q = 1:numel(ok)
        cln(q, :) = mode_class(bvn, real(V(:, ok(q))), ng, nd, nl, ma, 4 * nn, 1);
      end
      clas_nat{l} = cln;
    end
    cf.classification = clas;
    cf.classification_natural = clas_nat;
  end

  r = struct();
  r.cfsm = cf;
  r.load_factors_dense = lf_dense;
  r.name = cs.name;
  if isfield(cs, 'template'), r.template = cs.template; end
  r.E = E;
  r.nu = nu;
  r.actions = ac;
  r.constraints = constraints;
  r.springs = springs;
  r.neigs = neigs;
  r.node = node;
  r.elem = elem(:, 1:4);
  [Py, Mxx_y, Mzz_y, M11_y, M22_y] = yieldMP(node, 345, A, xcg, zcg, Ixx, Izz, Ixz, thetap, I11, I22, ac.unsymm);
  r.yield = struct('fy', 345, 'Py', Py, 'Mxx', Mxx_y, 'Mzz', Mzz_y, 'M11', M11_y, 'M22', M22_y);
  r.props = struct('A', A, 'xcg', xcg, 'zcg', zcg, 'Ixx', Ixx, 'Izz', Izz, 'Ixz', Ixz, 'thetap', thetap, 'I11', I11, 'I22', I22);
  r.bc = BC;
  r.lengths = lengths;
  r.m_all = m_all;
  r.strip1 = struct('k_local', full(kl), 'kg_local', full(kgl), 'k_global', full(kgk), 'kg_global', full(kgkg));
  r.K = K;
  r.Kg = Kg;
  r.load_factors = lf;
  r.mode1 = mode1;
  out{c} = r;
  printf('%s: %d lengths, %d nodes, lowest load factor %.6g\n', cs.name, numel(lengths), nn, min(cellfun(@(v) min([v Inf]), lf)));
end
fid = fopen(out_file, 'w');
fputs(fid, jsonencode(out));
fclose(fid);
