function out = test_helpers(action, varargin)
% TEST_HELPERS  Lightweight assertion utilities for the eVTOL test suite.
%
% Usage:
%   tr = test_helpers('init', 'test_name');
%   tr = test_helpers('assert', tr, condition, label, [diag_msg]);
%   tr = test_helpers('assert_near', tr, actual, expected, tol, label);
%   tr = test_helpers('assert_lt', tr, x, threshold, label);
%   test_helpers('report', tr);
%   tr.passed (bool), tr.npass, tr.ntotal

switch action
    case 'init'
        out.name   = varargin{1};
        out.npass  = 0;
        out.ntotal = 0;
        out.fails  = {};
        out.passed = true;

    case 'assert'
        tr = varargin{1};
        cond = varargin{2};
        label = varargin{3};
        diag = '';
        if numel(varargin) >= 4, diag = varargin{4}; end
        tr.ntotal = tr.ntotal + 1;
        if cond
            tr.npass = tr.npass + 1;
        else
            tr.passed = false;
            tr.fails{end+1} = sprintf('  [FAIL] %s%s', label, ...
                ifelse(isempty(diag), '', sprintf(' — %s', diag)));
        end
        out = tr;

    case 'assert_near'
        tr = varargin{1};
        actual = varargin{2};
        expected = varargin{3};
        tol = varargin{4};
        label = varargin{5};
        err = norm(actual - expected, 'inf');
        cond = err < tol;
        diag = sprintf('err=%.3e tol=%.3e', err, tol);
        out = test_helpers('assert', tr, cond, label, diag);

    case 'assert_lt'
        tr = varargin{1};
        x = varargin{2};
        threshold = varargin{3};
        label = varargin{4};
        cond = all(x(:) < threshold);
        diag = sprintf('value=%.3e thresh=%.3e', max(x(:)), threshold);
        out = test_helpers('assert', tr, cond, label, diag);

    case 'report'
        tr = varargin{1};
        if tr.passed
            tag = 'PASS';
        else
            tag = 'FAIL';
        end
        fprintf('  %-30s %s (%d/%d)\n', tr.name, tag, tr.npass, tr.ntotal);
        for k = 1:numel(tr.fails)
            fprintf('%s\n', tr.fails{k});
        end
        out = tr;

    otherwise
        error('test_helpers: unknown action %s', action);
end
end

function r = ifelse(c, a, b)
    if c, r = a; else, r = b; end
end
