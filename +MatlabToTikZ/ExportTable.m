% NOTE: if you want to use formulas in the headers or captions, you should
% use the double slash, i.e., write something like $\\phi_{m = 1}^{\\ast}$
%
function ExportTable(	strFileName,		...
						aacMatrix,			... either matrix of floats or cells (so to support also text)
						astrRowsHeaders,	... optional
						astrColsHeaders,	... optional
						iPrecision,			... optional
						strLaTeXLabel,		... optional
						strLaTeXCaption		) % optional
	%
	switch( nargin )
		%
		case 2
			astrRowsHeaders	= [];
			astrColsHeaders	= [];
			iPrecision		= 2;
			strLaTeXLabel	= [];
			strLaTeXCaption	= [];
		%
		case 3
			astrColsHeaders	= [];
			iPrecision		= 2;
			strLaTeXLabel	= [];
			strLaTeXCaption	= [];
		%
		case 4
			iPrecision		= 2;
			strLaTeXLabel	= [];
			strLaTeXCaption	= [];
		%
		case 5
			strLaTeXLabel	= [];
			strLaTeXCaption	= [];
		%
		case 6
			strLaTeXCaption	= [];
		%
	end;% switch on the number of arguments	
	%
	% in case somebody makes a mistake and passes a matrix of floats then convert it
	if( isfloat( aacMatrix ) )
		aacMatrix = num2cell( aacMatrix );
	end;%
	%
	% for readability
	iNumRows = numel( aacMatrix(:,1) );
	iNumCols = numel( aacMatrix(1,:) );
	%
	% check that there are no problems of dimensions
	if( numel(astrRowsHeaders) ~= 0 && numel(astrRowsHeaders) ~= iNumRows )
		warning( 'wrong number of rows headers -- I am going to ignore them' );
		astrRowsHeaders = [];
	end;%
	if( numel(astrColsHeaders) ~= 0 && numel(astrColsHeaders) ~= iNumCols )
		warning( 'wrong number of columns headers -- I am going to ignore them' );
		astrColsHeaders = [];
	end;%
	%
	% open the file
	fid = fopen(strFileName, 'w');
	%
	% write the preamble of the table
	fprintf(fid, '\\begin{table}[!htbp]\n');
	fprintf(fid, '\\centering\n');
	fprintf(fid, '\\begin{tabular}{');
	if( numel(astrRowsHeaders) ~= 0 )
		fprintf(fid, 'l|'); % put the space for the rows headers, in case
	end;%
	for iCol = 1:iNumCols;
		fprintf(fid, 'c');
	end;%
	fprintf(fid, '}\n');
	%
	if( numel(astrColsHeaders) ~= 0 )
		if( numel(astrRowsHeaders) ~= 0 )
			fprintf(fid, '& '); % put the space for the rows headers, in case
		end;%
		for iCol = 1:iNumCols-1
			fprintf(fid, '%s \t & ', astrColsHeaders{iCol} );
		end;%
		fprintf(fid, '%s \t \\\\ \n', astrColsHeaders{end} );
		fprintf(fid, '\\hline\n');
	end;%
	%
	% write the actual table
	for iRow = 1:iNumRows;
		%
		if( numel(astrRowsHeaders) ~= 0 )
			fprintf(fid, astrRowsHeaders{iRow});
			fprintf(fid, '\t & ');
		end;%
		%
		for iCol = 1:iNumCols-1
			if( isstring(aacMatrix{iRow, iCol}) || ischar(aacMatrix{iRow, iCol}) )
				fprintf(fid, '%s \t & ', aacMatrix{iRow, iCol} );
			elseif( isfloat(aacMatrix{iRow, iCol}) )
				fprintf(fid, '%.*f \t & ', iPrecision, aacMatrix{iRow, iCol} );
			elseif( isinteger(aacMatrix{iRow, iCol}) )
				fprintf(fid, '%d \t & ', aacMatrix{iRow, iCol} );
			else%
				error('unsupported data type');
			end;%
		end;%
		%
		if( isstring(aacMatrix{iRow, end}) || ischar(aacMatrix{iRow, end}) )
			fprintf(fid, '%s \t \\\\ \n', aacMatrix{iRow, end} );
		elseif( isfloat(aacMatrix{iRow, end}) )
			fprintf(fid, '%.*f \t \\\\ \n', iPrecision, aacMatrix{iRow, end} );
		elseif( isinteger(aacMatrix{iRow, end}) )
			fprintf(fid, '%d \t \\\\ \n', aacMatrix{iRow, end} );
		else%
			error('unsupported data type');
		end;%
		%
	end;%
	%
	% write the ending of the tabular environment
	fprintf(fid, '\\end{tabular}\n');
	%
	% if there is a caption, then write it down
	if( ~strcmp(strLaTeXCaption, []) )
		fprintf(fid, '\\caption{');
		fprintf(fid, sprintf('%s', strLaTeXCaption));
		fprintf(fid, '}\n');
	end;%
	%
	% if there is a label, then write it down
	if( ~strcmp(strLaTeXLabel, []) )
		fprintf(fid, '\\label{');
		fprintf(fid, sprintf('%s', strLaTeXLabel));
		fprintf(fid, '}\n');
	end;%
	%
	% write the ending of the table environment
	fprintf(fid, '\\end{table}\n');
	%
	% close the file
	fclose(fid);
	clear fid;
	fprintf('Table %s exported.\n', strFileName);
	%
end % function

