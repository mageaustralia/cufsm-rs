function varargout = eigs(varargin)
% Oracle shim: stripmain.m calls eigs(K, Kg, N, 'SM', 'Display', 0), the MATLAB name-value form.
% Octave's eigs takes an options struct instead, so the pair is dropped (display is off by
% default) and the call passed to Octave's own eigs, kept here as eigs_octave.m.
args = varargin;
i = find(cellfun(@(x) ischar(x) && strcmpi(x, 'Display'), args), 1);
if ~isempty(i)
  args(i:i+1) = [];
end
[varargout{1:max(nargout, 1)}] = eigs_octave(args{:});
end
