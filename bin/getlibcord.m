function libcord = getlibcord(fileloc)

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 5, "Encoding", "UTF-8");

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["cord384", "ORF", "Gene", "Var4", "Var5"];
opts.SelectedVariableNames = ["cord384", "ORF", "Gene"];
opts.VariableTypes = ["string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["cord384", "ORF", "Gene", "Var4", "Var5"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["cord384", "ORF", "Gene", "Var4", "Var5"], "EmptyFieldRule", "auto");

% Import the data
libcord = readtable(fileloc, opts);


%% Clear temporary variables
clear opts
end