% Extracts the models and saved results from CUFSM's AISI Direct Strength Method Design Guide
% (2006) example files, as JSON for the Rust tests:
%
%   octave-cli --no-gui extract_dsm_guide.m <dir of .mat files> <out.json>
%
% Each file holds the model (prop, node, elem), the half-wavelengths, and the load factors CUFSM
% computed then (curve: lengths x [length lf] x modes, simply supported, one term). Files with
% springs (the pre-4.3 format), no results, or constraints of another shape are skipped.
args = argv();
d = args{1};
out_file = args{2};
files = dir(fullfile(d, '*.mat'));
out = {};
for f = 1:numel(files)
  s = load(fullfile(d, files(f).name));
  if ~all(isfield(s, {'prop', 'node', 'elem', 'lengths', 'curve'})) || isempty(s.curve) || ~isnumeric(s.curve)
    printf('skip %s: no model or results\n', files(f).name); continue;
  end
  if isfield(s, 'springs') && ~isempty(s.springs) && any(s.springs(:) ~= 0)
    printf('skip %s: springs\n', files(f).name); continue;
  end
  cons = [];
  if isfield(s, 'constraints') && ~isempty(s.constraints) && any(s.constraints(:) ~= 0)
    cons = s.constraints;
    if size(cons, 2) < 5, printf('skip %s: constraints\n', files(f).name); continue; end
  end
  c = s.curve;
  if ndims(c) ~= 3 && ~(ismatrix(c) && size(c, 2) == 2)
    printf('skip %s: curve shape\n', files(f).name); continue;
  end
  lf = squeeze(c(:, 2, :));
  if isvector(lf), lf = lf(:); end
  r = struct();
  r.name = files(f).name;
  r.prop = s.prop;
  r.node = s.node;
  r.elem = s.elem;
  r.lengths = s.lengths(:)';
  r.curve_lengths = c(:, 1, 1)';
  r.load_factors = lf;
  r.constraints = cons;
  % cFSM spaces the run was restricted to, if any: [glob dist local other] counts switched on.
  g = [0 0 0 0];
  if isfield(s, 'GBTcon') && isstruct(s.GBTcon) && isfield(s.GBTcon, 'glob')
    g = [sum(s.GBTcon.glob) sum(s.GBTcon.dist) sum(s.GBTcon.local) sum(s.GBTcon.other)];
    r.gbt_sizes = [numel(s.GBTcon.glob) numel(s.GBTcon.dist) numel(s.GBTcon.local) numel(s.GBTcon.other)];
  end
  r.gbt = g;
  out{end + 1} = r;
  printf('%s: %d nodes, %d lengths, %d modes\n', files(f).name, size(s.node, 1), size(lf, 1), size(lf, 2));
end
fid = fopen(out_file, 'w');
fputs(fid, jsonencode(out));
fclose(fid);
