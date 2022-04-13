

function SW_fun(input_file, similarity_matrix, openGap, extGap)

%   This function will align two sequences based on Smith-Waterman
%   algorithm to create the best alignment possible. The algorithm will use
%   an afine gap penalty
%
%   The result of this function will be a txt file with the scoring matrix
%   and best alignment.
%
%   Usage: 
%           SW_fun(input_file, similarity_matrix, openGap, extGap)
%           SW_fun('input.txt', 'blosum62.txt', 2, 1)
%           SW_fun('input.txt', 'blosum62.txt')
%
%   input_file -->  should be in txt format. If the input file is not in the
%                   working directory then the complete directory location 
%                   is inputed as the input file. The input file should be
%                   a character when used in the function
%
%   similarity_matrix -->  should be in txt format. Similar to input file
%                          if the similarity matrix is not in the working
%                          directory then the complete location is
%                          provided. The name should be a character when
%                          used in the function
%
%   openGap  --> is a number digit. This will specify the opening gap cost
%                of an alignment function. if not specified then the
%                default value of 2 will be used
%
%   extGap  --> is a number digit. This will specify the extension gap cost
%                of an alignment function. if not specified then the
%                default value of 1 will be used
%
%
%   Author: Alaaddin Ibrahimy
%   Date:   04/13/2022

%   MATLAB versions 2019a or later
%   export PATH=/Applications/MATLAB_R2021b.app/bin/:$PATH
%--------------------------------------------------------------------------


%%  Input files

inputname = input_file;
simmat = similarity_matrix;

%--------------------------------------------------------------------------
% Affine Gap equation --> Wk = uk + v  (where, u = extGap, v = openGap)
% Gap values should be positive because we subtract them later in the
% equation
if nargin < 3
    openGap = 2;                   % Opening gap value
    extGap = 1;                    % Extension gap value
end

%----- Making sure that the gap penalties are positive numbers
if openGap < 0
    openGap = - openGap;
end
if extGap < 0
    extGap = -extGap;
end


%%  Reading the similarity and input file with sequences
%--------------------------------------------------------------------------
simdata = readcell(simmat);      % Similarity matrix
simdata(1,2:end) = simdata(1,1:end-1);
simdata{1,1} = [];

%---- Loading input sequences from text file
inputdata = readcell(inputname); % input sequence data
seq_1 = inputdata{1};            % sequence 1 to be shown in Top row
n = length(seq_1);               % Size of the sequence 1
seq_2 = inputdata{2};            % Sequence 2 to be shown in the right col
m = length(seq_2);               % Size of the sequence 2

%%  Creating scoring matrix
%--------------------------------------------------------------------------

H = zeros(m+1, n+1);    % initializing scoring matrix

flag = cell(m+1, n+1);  % Initializing a flag matrix to define the direction
for i = 1 :size(flag,1)
    for j = 1:size(flag,2)
        flag{i,j} = 'nan';
    end
end


for i = 2:m+1          % looping in row (seq_2)
    for j = 2:n+1      % looping in coloum (seq_1)
        %------------------------------------------------------------------
        %   Calculating the diagonal value
        %   First find the row and col value corresponding to the
        %   simmulaiton matrix and then add to the diagnoal value
        temSeq_1 = seq_1(j-1);      % Sequence 1 letter
        temSeq_2 = seq_2(i-1);      % Sequence 2 letter
        row = find([simdata{:,1}] == temSeq_1) +1;     % row number that has letter for seq_1 in similarity matrix
        col = find([simdata{1,:}] == temSeq_2) +1;     % Col number that has letter for seq_2 in similarity matrix
        match = H(i-1,j-1) + simdata{row,col};         % Diagonal value + similarity matrix score

        %------------------------------------------------------------------
        %   Calculating deletion and insertion matrix
        %   looping over the rows and column values to find the max minus
        %   the gap penalty
        for k = 1:i-1
            deletion(k) = H(i-k,j) - (openGap + extGap * (k-1));
        end
        deletion = max(deletion);   % Finding the max of deletion
        for k = 1:j-1
            insertion(k) = H(i,j-k) - (openGap + extGap * (k-1));
        end
        insertion = max(insertion); % Finding the max of insertion

        %------------------------------------------------------------------
        %   Assigning the Max value based on the Diag, Up, Left, Zero to
        %   the cell we are looping over
        H(i,j) = max([match, deletion, insertion, 0]);

        %------------------------------------------------------------------
        %   Filling out the flag matrix to see which way it was selected
        tempmax = max([0, match, deletion, insertion]);
        if tempmax == 0             % if zero was the max
            flag{i,j} = 'nan';
        elseif tempmax == match     % if diagonal was the max
            flag{i,j} = 'Diag';
        elseif tempmax == deletion  % if Up was the max
            flag{i,j} = 'Up';
        elseif tempmax == insertion % if left was the max
            flag{i,j} = 'Left';
        end

    end
end

Score = max(max(H));        % Highest score value

%%  Back Tracing the scoring matrix
%--------------------------------------------------------------------------

[row, col] = find(H == max(max(H)));    % The location of the maximum score

seq_1_align = [];    % initializing aligned sequence 1
seq_2_align = [];    % initializing aligned sequence 2
lines = [];          % initializing the lines for the location of aligned sequences

% creating row and column index to traceback
r = row;
c = col;
while ~contains(flag{r,c}, 'nan')
    direction = flag{r,c};
    if contains(direction, 'Diag')        % if traceback was from diagonal

        seq_1_align = [seq_1(c-1), seq_1_align];
        seq_2_align = [seq_2(r-1),seq_2_align];

        if seq_1(c-1) == seq_2(r-1)
            lines = ['|', lines];
        else
            lines = [' ', lines];
        end
        r = r-1;
        c = c-1;

    elseif contains(direction, 'Left')   % if traceback was from Left

        seq_1_align = [seq_1(c-1), seq_1_align];
        seq_2_align = ['-', seq_2_align];
        lines = [' ', lines];
        c = c-1;

    elseif contains(direction, 'Up')     % if traceback was from Up
        seq_1_align = ['-', seq_1_align];
        seq_2_align = [seq_2(r-1), seq_2_align];
        lines = [' ', lines];
        r = r-1;
    end
end

%%  Padding the beginning and end of sequences
%------- Adding the begining part of each sequence
seq_1_align = [seq_1(1:c-1),'(', seq_1_align, ')', seq_1(col:end)];     % Begining and end part of aligned sequence 1
seq_2_align = [seq_2(1:r-1), '(', seq_2_align,')', seq_2(row:end)];     % Begining and end part of aligned sequence 2
lines = [' ',lines, ' '];

s1lenpad = length(seq_1(col:end));    % length of the padding letter for sequence 1
s2lenpad = length(seq_2(row:end));    % length of the padding letter for sequence 2
s1lenbeg = length(seq_1(1:c-1));      % length of the begining letter for sequence 1
s2lenbeg = length(seq_2(1:r-1));      % length of the begining letter for sequence 2

%------ Adding the spaces to begining and end of the lines
for i = 1: max([s1lenpad; s2lenpad])
    lines = [lines, ' '];
end
for i = 1 : max([s1lenbeg; s2lenbeg])
    lines = [' ', lines];
end


%------- Adding space at the end of the sequences
if s1lenpad > s2lenpad  % if sequence 1 has padding, add space for sequence 2
    for i = 1 : (s1lenpad - s2lenpad)
        seq_2_align = [seq_2_align, ' '];
    end
else         % if sequence 2 has padding, add space to sequence 1
    for i = 1 : (s2lenpad - s1lenpad)
        seq_1_align = [seq_1_align, ' '];
    end
end

%------- Adding space at the beginning of the sequences
if s1lenbeg > s2lenbeg  % if sequence 1 has padding, add space for sequence 2
    for i = 1 : (s1lenbeg - s2lenbeg)
        seq_2_align = [' ', seq_2_align];
    end
else         % if sequence 2 has padding, add space to sequence 1
    for i = 1 : (s2lenbeg - s1lenbeg)
        seq_1_align = [' ', seq_1_align];
    end

end


%--- Aligned sequnces in one character vector
x = [seq_1_align;...
    lines;...
    seq_2_align];


%%  Writting the output data into one file to export
%--------------------------------------------------------------------------
%--------   Titles
OutputFile{1,1} = '-----------';
OutputFile{2,1} = '|Sequences|';
OutputFile{3,1} = '-----------';
OutputFile{4,1} = 'sequence1';
OutputFile{5,1} = seq_1;
OutputFile{6,1} = 'sequence2';
OutputFile{7,1} = seq_2;
OutputFile{8,1} = '--------------';
OutputFile{9,1} = '|Score Matrix|';
OutputFile{10,1} ='--------------';

%----- Score Matrix
for i = 1 : size(H,1)       % Rows of Score matrix H
    for j = 1 : size(H,2)   % Columns of Score matrix H
        OutputFile{i+11,j+1} = H(i,j);
    end
end

%------- Sequence 1 titles
for j = 3 : length(seq_1) + 2
    OutputFile{11,j} = seq_1(j-2);
end

%------- Sequence 2 titles
for i = 13 : length(seq_2) + 12
    OutputFile{i,1} = seq_2(i-12);
end

%--------- Alignment and Score
mm = length(seq_2) + 12;
OutputFile{mm+1,1} = '----------------------';
OutputFile{mm+2,1} = '|Best Local Alignment|';
OutputFile{mm+3,1} = '----------------------';
OutputFile{mm+4,1} = ['Alignment Score:',num2str(Score)];
OutputFile{mm+5,1} = 'Alignment Results:';
OutputFile{mm+6,1} = seq_1_align;
OutputFile{mm+7,1} = lines;
OutputFile{mm+8,1} = seq_2_align;


%---- Writing the file
writecell(OutputFile,'output.txt', 'delimiter','tab')









