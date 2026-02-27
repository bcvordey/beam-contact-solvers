function latexRHS = psiToLatexRHS(psiFun)
% Convert a function handle (e.g., @(x) -0.02 + 0.01*sin(12*pi*x).^2)
% to a LaTeX RHS string (e.g., -0.02 + 0.01\sin^2(12\pi x)).
% Prints 's' if verbose==true.
% escapeForSprintf (default true) doubles backslashes for sprintf/title.

    s = func2str(psiFun);                    % '@(x)-0.02+0.01*sin(12*pi*x).^2'
    s = regexprep(s,'^\s*@\(\w+\)\s*','');   % strip '@(x)'

    % elementwise -> mathy
    s = strrep(s,'.^','^');
    s = strrep(s,'./','/');
    s = strrep(s,'.*','*');

    % sin/cos/... squared: sin(...) .^ 2 -> \sin^2(...)
    s = regexprep(s,'\b(sin|cos|tan|sinh|cosh|tanh)\(\s*([^)]+?)\s*\)\s*\.\^\s*2','\\$1^2($2)');

    % functions -> LaTeX
    s = regexprep(s,'\b(sin|cos|tan|sinh|cosh|tanh|exp)\(','\\$1(');
    s = regexprep(s,'\bsqrt\(\s*([^)]+)\s*\)','\\sqrt{$1}');

    % PI via star patterns: *pi*, *pi, pi*
    s = regexprep(s,'\s*\*\s*pi\s*\*\s*','\\pi ');
    s = regexprep(s,'\s*\*\s*pi\b','\\pi');
    s = regexprep(s,'\bpi\s*\*\s*','\\pi ');

    % remove remaining '*'
    s = regexprep(s,'\s*\*\s*',' ');

    % tidy
    s = regexprep(s,'\\pi(?=[A-Za-z])','\\pi ');
    s = regexprep(s,'\s+',' ');
    s = strtrim(s);
    s = strrep(s,'\','\\'); 

    fprintf('[psiToLatexRHS] pre-escape: %s\n', s);
 
    latexRHS = s;
end
