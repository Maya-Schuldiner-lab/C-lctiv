function SGDDB = SGD_DB_import(filename)

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = "|";

% Specify column names and types
opts.VariableNames = ["ORF", "Gene", "description", "protein_sequence", "DNA_sequence"];
opts.VariableTypes = ["string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["ORF", "Gene", "description", "protein_sequence", "DNA_sequence"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["ORF", "Gene", "description", "protein_sequence", "DNA_sequence"], "EmptyFieldRule", "auto");

% Import the data
SGDDB = readtable(filename, opts);

end