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

  r = struct();
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
