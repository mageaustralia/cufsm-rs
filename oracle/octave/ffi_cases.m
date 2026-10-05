% The reference values for the C interface's console exports (tests/ffi_exports.rs) that the main
% fixture does not hold, computed with CUFSM's own MATLAB code on the main fixture's models.
%
%   octave-cli --no-gui --norc ffi_cases.m <cufsm-root> <cufsm_octave.json> <out.json>
%
% (under micromamba as run_octave.sh does). <cufsm-root> is a checkout of CUFSM
% (https://github.com/thinwalled/cufsm-git); nothing in it is modified.
%
% 1. signature_ss.m's 100 half-wavelengths for every model. signature_ss.m is run as it is, with
%    a stand-in stripmain.m (in a temporary folder) that returns the lengths it is given instead
%    of solving at them: only the lengths are wanted here.
% 2. Load factors with a different, non-contiguous set of longitudinal terms at each length, on
%    S-S, C-C and C-F models, one with springs: CUFSM's assembly and constraint basis solved with
%    eig(), as run_cases.m's stage 3b.
% 3. Load factors restricted to unions of cFSM spaces (G+D, D+L, L+O, G+D+L), as run_cases.m's
%    stage 4 restricts to one.
% 4. The classification of the first three modes under every orth and norm option CUFSM offers
%    (orth 1, 2, 3 with norm 0, 1, 2, 3; the ST O space), as run_cases.m's stage 4 does for some.

args = argv();
root = args{1};
in_file = args{2};
out_file = args{3};
here = fileparts(mfilename('fullpath'));

cases = jsondecode(fileread(in_file));
if ~iscell(cases), cases = num2cell(cases); end
[~, sha] = system(['git -C "' root '" rev-parse --short HEAD']);
res = struct();
res.about = ['CUFSM ' strtrim(sha) ' under Octave ' version() ', from oracle/octave/ffi_cases.m on the models of cufsm_octave.json'];
res.cufsm_commit = strtrim(sha);

% Stage 1: signature_ss.m's lengths, with the stand-in stripmain.m first on the path.
stub = tempname(); mkdir(stub);
fid = fopen(fullfile(stub, 'stripmain.m'), 'w');
fprintf(fid, '%s\n', ...
  'function [curve, shapes] = stripmain(prop, node, elem, lengths, springs, constraints, GBTcon, BC, m_all, neigs)', ...
  'curve = cell(numel(lengths), 1);', ...
  'for i = 1:numel(lengths), curve{i} = [lengths(i) 0]; end', ...
  'shapes = {};', ...
  'end');
fclose(fid);
addpath(fullfile(root, 'helpers'));
addpath(stub);
sig = cell(numel(cases), 1);
for c = 1:numel(cases)
  cs = cases{c};
  prop = [100 cs.E cs.E cs.nu cs.nu cs.E / (2 * (1 + cs.nu))];
  elem = cs.elem; elem(:, 5) = 100;
  curve = signature_ss(prop, cs.node, elem, struct('glob', 0, 'dist', 0, 'local', 0, 'other', 0));
  sig{c} = struct('name', cs.name, 'lengths', cellfun(@(r) r(1), curve)');
end
res.signature_lengths = sig;
rmpath(stub);

addpath(fullfile(root, 'analysis'));
addpath(fullfile(root, 'analysis', 'cFSM'));

% Assembled K and Kg at length al with terms ma, springs added as stripmain.m adds them.
function [KL, KgL] = assemble_all(node, elem, E, nu, al, BC, ma, springs)
  nn = size(node, 1);
  elprop = elemprop(node, elem, nn, size(elem, 1));
  KL = sparse(zeros(4 * nn * numel(ma))); KgL = KL;
  for i = 1:size(elem, 1)
    ti = elem(i, 4); bi = elprop(i, 2);
    k_l = klocal(E, E, nu, nu, E / (2 * (1 + nu)), ti, al, bi, BC, ma);
    kg_l = kglocal(al, bi, node(elem(i, 2), 8) * ti, node(elem(i, 3), 8) * ti, BC, ma);
    [k, kg] = trans(elprop(i, 3), k_l, kg_l, ma);
    [KL, KgL] = assemble(KL, KgL, k, kg, elem(i, 2), elem(i, 3), nn, ma);
  end
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
end

% The positive load factors of Kff phi = lambda Kgff phi, smallest first, at most neigs.
function lf = positive(Kff, Kgff, neigs)
  ev = eig(full(Kff), full((Kgff + Kgff') / 2));
  ev = real(ev(abs(imag(ev)) < 1e-5 * abs(real(ev)) & real(ev) > 0));
  ev = sort(ev);
  lf = ev(1:min(neigs, numel(ev)))';
end

function [springs, constraints] = extras(cs)
  springs = 0; constraints = 0;
  if isfield(cs, 'springs') && ~isempty(cs.springs) && ~isequal(cs.springs, 0)
    springs = cs.springs; if isvector(springs), springs = springs(:)'; end
  end
  if isfield(cs, 'constraints') && ~isempty(cs.constraints) && ~isequal(cs.constraints, 0)
    constraints = cs.constraints; if isvector(constraints), constraints = constraints(:)'; end
  end
end

function c = by_name(cases, name)
  for k = 1:numel(cases)
    if strcmp(cases{k}.name, name), c = cases{k}; return; end
  end
  error('no case %s', name);
end

% Stage 2: a different, non-contiguous set of terms at each length.
picks = {'lipped-c S-S compression', 'S-S', [300 1200 4000], {[1 3 5], [2 4 6 8], [1 2 3 7]}; ...
         'lipped-c C-C compression', 'C-C', [300 1200 4000], {[1 3 5], [2 4 6 8], [1 2 3 7]}; ...
         'lipped-c C-C springs', 'C-C', [600 2500], {[1 2 4], [3 5 7 9]}; ...
         'lipped-c C-F compression', 'C-F', [300 1200 4000], {[1 3], [1 2 5 6], [2 3 4]}};
terms = cell(size(picks, 1), 1);
for p = 1:size(picks, 1)
  cs = by_name(cases, picks{p, 1});
  [springs, constraints] = extras(cs);
  node = cs.node; elem = cs.elem; nn = size(node, 1);
  BC = picks{p, 2}; lengths = picks{p, 3}; m_all = picks{p, 4};
  lfs = cell(numel(lengths), 1);
  for l = 1:numel(lengths)
    ma = msort({m_all{l}}){1};
    [KL, KgL] = assemble_all(node, elem, cs.E, cs.nu, lengths(l), BC, ma, springs);
    if constr_BCFlag(node, constraints) == 0
      R = speye(4 * nn * numel(ma));
    else
      R = null(null(constr_user(node, constraints, ma)')');
    end
    lfs{l} = positive(R' * KL * R, R' * KgL * R, 5);
  end
  terms{p} = struct('name', cs.name, 'bc', BC, 'lengths', lengths, 'm_all', {m_all}, 'load_factors', {lfs});
end
res.terms = terms;

% Stages 3 and 4 on the cFSM cases whose spaces are all non-empty and whose bases are unique
% enough to compare (the plain channel has no D space; the doubly symmetric I-section's axial
% basis is not unique).
cfsm_names = {'cfsm lipped-c', 'cfsm lipped-c rounded', 'cfsm lipped-z', 'cfsm hat', 'cfsm lipped-c fixed and sprung'};
unions = {'GD', [1 1 0 0]; 'DL', [0 1 1 0]; 'LO', [0 0 1 1]; 'GDL', [1 1 1 0]};
opts = [1 0; 1 1; 1 2; 1 3; 2 0; 2 1; 2 2; 2 3; 3 0; 3 1; 3 2; 3 3];
cf = cell(numel(cfsm_names), 1);
for c = 1:numel(cfsm_names)
  cs = by_name(cases, cfsm_names{c});
  [springs, constraints] = extras(cs);
  node = cs.node; elem = cs.elem; elem(:, 5) = 100; nn = size(node, 1);
  prop = [100 cs.E cs.E cs.nu cs.nu cs.E / (2 * (1 + cs.nu))];
  lengths = cs.lengths(:)'; BC = cs.bc; ma = 1;
  rec = struct('name', cs.name);
  for u = 1:size(unions, 1)
    fl = unions{u, 2};
    lfs = cell(numel(lengths), 1);
    for l = 1:numel(lengths)
      al = lengths(l);
      [bv, ng, nd, nl] = base_column(node, elem, prop, al, BC, ma);
      R = mode_select(bv, ng, nd, nl, ones(1, ng) * fl(1), ones(1, nd) * fl(2), ones(1, nl) * fl(3), ones(1, 4 * nn - ng - nd - nl) * fl(4), 4 * nn, ma);
      [KL, KgL] = assemble_all(node, elem, cs.E, cs.nu, al, BC, ma, springs);
      if constr_BCFlag(node, constraints) ~= 0
        R = null([null(R') null(constr_user(node, constraints, ma)')]');
      end
      if isempty(R) || size(R, 2) == 0, lfs{l} = []; continue; end
      lfs{l} = positive(R' * KL * R, R' * KgL * R, 5);
    end
    rec.(['lf_' unions{u, 1}]) = lfs;
  end
  % Classification: the bare strips' modes (no springs or constraints, as run_cases.m's stage 4).
  if isequal(springs, 0) && isequal(constraints, 0) && all(all(node(:, 4:7) == 1))
    for o = 1:size(opts, 1)
      orth = opts(o, 1); nrm = opts(o, 2);
      clas = cell(numel(lengths), 1);
      for l = 1:numel(lengths)
        al = lengths(l);
        [KL, KgL] = assemble_all(node, elem, cs.E, cs.nu, al, BC, ma, 0);
        [V, D] = eig(full(KL), full((KgL + KgL') / 2));
        ev = diag(D); ok = find(abs(imag(ev)) < 1e-5 * abs(real(ev)) & real(ev) > 0);
        [~, ord] = sort(real(ev(ok))); ok = ok(ord); ok = ok(1:min(3, numel(ok)));
        [bv, ng, nd, nl] = base_column(node, elem, prop, al, BC, ma);
        bvu = base_update(1, nrm, bv, al, ma, node, elem, prop, ng, nd, nl, BC, 1, orth);
        clq = zeros(numel(ok), 4);
        for q = 1:numel(ok)
          clq(q, :) = mode_class(bvu, real(V(:, ok(q))), ng, nd, nl, ma, 4 * nn, 1);
        end
        clas{l} = clq;
      end
      rec.(sprintf('classification_orth%d_norm%d', orth, nrm)) = clas;
    end
  end
  cf{c} = rec;
  printf('%s done\n', cs.name);
end
res.cfsm = cf;

fid = fopen(out_file, 'w');
fputs(fid, jsonencode(res));
fclose(fid);
