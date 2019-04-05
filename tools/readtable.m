function t = readtable(filename,varargin)
%READTABLE Create a table by reading from a file.
%   Use the READTABLE function to create a table by reading column-oriented data
%   from a file. READTABLE can automatically determine the file format from its
%   extension as described below.
%
%   T = READTABLE(FILENAME) creates a table by reading from the file FILENAME,
%   and determines the file format from its extension. The extension must be
%   one of those listed below.
%
%   T = READTABLE(FILENAME,'FileType',FILETYPE) specifies the file type, where
%   FILETYPE is one of 'text' or 'spreadsheet'.
%
%   T = READTABLE(FILENAME,OPTS) creates a table from the file FILENAME
%   using the supplied ImportOptions OPTS. OPTS specifies variable names,
%   selected variable names, variable types, and other information
%   regarding the location of the data.
%
%   For example, import a sub-set of the data in a file:
%       opts = detectImportOptions('patients.xls')
%       opts.SelectedVariableNames = {'Systolic','Diastolic'}
%       T = readtable('patients.xls',opts)
%
%   READTABLE reads data from different file types as follows:
%
%   .txt, .dat, .csv:  Delimited text file.
%
%          Reading from a delimited text file creates one variable in T for
%          each column in the file. Variable names are taken from the first row
%          of the file. By default, the variables created are either double,
%          if the entire column is numeric, or cell array of strings, if any
%          element in a column is not numeric. READTABLE converts empty fields
%          in the file to either NaN (for a numeric variable) or the empty string
%          (for a string-valued variable). Insignificant whitespace in the file
%          is ignored.
%
%          Use the following optional parameter name/value pairs to control how
%          data are read from a delimited text file:
%
%          'Delimiter'     The delimiter used in the file. Can be any of ' ',
%                          '\t', ',', ';', '|' or their corresponding string
%                          names 'space', 'tab', 'comma', 'semi', or 'bar'.
%                          If unspecified, READTABLE detects the delimiter 
%                          automatically.
%
%          'ReadVariableNames'  A logical value that specifies whether or not the
%                          first row (after skipping HeaderLines) of the file is
%                          treated as variable names. If unspecified, READTABLE 
%                          detects the presence of variable names automatically.
%
%          'ReadRowNames'  A logical value that specifies whether or not the
%                          first column of the file is treated as row names.
%                          Default is false. If the 'ReadVariableNames' and
%                          'ReadRowNames' parameter values are both true, the
%                          name in the first column of the first row is saved
%                          as the first dimension name for the table.
%
%          'TreatAsEmpty'  One or more strings to be treated as the empty string
%                          in a numeric column. This may be a character string,
%                          or a cell array of strings. Table elements
%                          corresponding to these are set to NaN. 'TreatAsEmpty'
%                          only applies to numeric columns in the file, and
%                          numeric literals such as '-99' are not accepted.
%
%          'HeaderLines'   The number of lines to skip at the beginning of the
%                          file. If unspecified, READTABLE detects the number of 
%                          lines to skip automatically.
%
%          'Format'        A format string to define the columns in the file, as
%                          accepted by the TEXTSCAN function. If you specify 'Format',
%                          you may also specify any of the parameter name/value pairs
%                          accepted by the TEXTSCAN function. Type "help textscan" for
%                          information about format strings and additional parameters.
%                          Specifying the format can significantly improve speed for
%                          some large files. If unspecified, READTABLE detects the
%                          format automatically.
%
%          'DateLocale'    A string specifying the locale that READTABLE uses to interpret
%                          month and day names in datetime strings read with the %D format
%                          specifier. LOCALE must be a string in the form xx_YY. See the
%                          documentation for DATETIME for more information.
%
%          'FileEncoding'  A string that specifies the character 
%                          encoding scheme associated with the file. It must 
%                          be the empty string ('') or a name or alias for 
%                          an encoding scheme. READTABLE supports the same  
%                          character encodings as FOPEN. Some examples are  
%                          'UTF-8', 'latin1', 'US-ASCII', and 'Shift_JIS'. 
%                          If the 'FileEncoding' parameter value is unspecified 
%                          or is the empty string (''), then READTABLE uses   
%                          your system's default encoding scheme is used.
%
%                          See the documentation for FOPEN for more information.
%
%          'TextType'      The output type of text variables. Text variables are 
%                          those with %s, %q, or %[...] formats. It can
%                          have either of the following values:
%                             'char'   - Return text as a cell array of character vectors.
%                             'string' - Return text as a string array.
%
%          'DatetimeType'  The output type of datetime variables. The
%                          possible values are
%                              'datetime'     - Return date and time data as MATLAB datetimes.
%                              'text'         - Return date and time data as text.
%
%   .xls, .xlsx, .xlsb, .xlsm, .xltm, .xltx, .ods:  Spreadsheet file.
%
%          Reading from a spreadsheet file creates one variable in T for each column
%          in the file. By default, the variables created are either double, or cell
%          array of strings. Variable names are taken from the first row of the
%          spreadsheet.
%
%          Use the following optional parameter name/value pairs to control how
%          data are read from a spreadsheet file:
%
%          'ReadVariableNames'  A logical value that specifies whether or not the
%                          first row of the specified region of the file is treated
%                          as variable names. Default is true.
%
%          'ReadRowNames'  A logical value that specifies whether or not the first
%                          column of specified region of the file is treated as row
%                          names. Default is false. If the 'ReadVariableNames'
%                          and 'ReadRowNames' parameter values are both true, the
%                          name in the first column of the first row is saved as
%                          the first dimension name for the table.
%
%          'TreatAsEmpty'  One or more strings to be treated as an empty cell
%                          in a numeric column. This may be a character string,
%                          or a cell array of strings. Table elements
%                          corresponding to these are set to NaN. 'TreatAsEmpty'
%                          only applies to numeric columns in the file, and
%                          numeric literals such as '-99' are not accepted.
%
%          'Sheet'         The sheet to read, specified as a string that contains
%                          the worksheet name, or a positive integer indicating the
%                          worksheet index.
%
%          'Range'         A string that specifies a rectangular portion of the
%                          worksheet to read, using the Excel A1 reference
%                          style. If the spreadsheet contains figures or other
%                          non-tabular information, you should use the 'Range'
%                          parameter to read only the tabular data. By default,
%                          READTABLE reads data from a spreadsheet contiguously
%                          out to the right-most column that contains data,
%                          including any empty columns that precede it. If the
%                          spreadsheet contains one or more empty columns
%                          between columns of data, use the 'Range' parameter to
%                          specify a rectangular range of cells from which to
%                          read variable names and data. An empty string ('')
%                          signifies all data in the sheet.
%
%          'Basic'         A logical value specifying whether or not to read the
%                          spreadsheet in basic mode. Basic mode is the default
%                          for systems without Excel for Windows installed.
%                          Basic mode may read files where there may be live
%                          updates (e.g. formula evalauation or plugins) to the 
%                          data differently. If your files do not have this 
%                          requirement, basic mode can be used.
%
%          'TextType'      The output type of text variables. It can have
%                          either of the following values:
%                             'char'   - Return text as a cell array of character vectors.
%                             'string' - Return text as a string array.
%
%          'DatetimeType'  The output type of datetime variables. The
%                          possible values are
%                              'datetime'     - Return date and time data as MATLAB datetimes.
%                              'text'         - Return date and time data as text.
%                              'exceldatenum' - Return date and time data as Excel serial day date numbers.
%
%   Parameters which are also accepted with import options. These may have
%   slightly different behavior when used with import options:
%
%       T = readtable(FILENAME, OPTS, ...parameters...)
%
%         'ReadVariableNames' true - Reads the variable names from
%                           opts.VariableNamesRange/VariableNamesLine location.
%                           false - Uses variable names from the import options. 
%                                       
%         'ReadRowNames'    true - Reads the row names from the
%                           opts.RowNamesRange/RowNameColumn location. 
%                           false - Does not import row names.
%
%       Text only parameters:
%         'DateLocale'      Override the locale used when importing dates.
%         'FileEncoding'    Override the encoding defined in import options.
%
%       Spreadsheet only parameters:
%         'Sheet'            Override the sheet value in the import options.
%         'Basic'            No difference.
%
%
%   See also WRITETABLE, TABLE, DETECTIMPORTOPTIONS, TEXTSCAN.

%   Copyright 2012-2016 The MathWorks, Inc.

t = table.readFromFile(filename,varargin);
