% ExportGraph.m
%
%
% inputs:
%
% * strFileName          = full path
% * tGraph, containing
%   * aafNodesPositions  = matrix with each row being the x and y of the nodes
%   * aaiAdjacencyMatrix = adjacency matrix
% * strCommunicationKind = 'directed' | 'undirected'
% * aacNodesStatus       = **NOT REQUIRED**, if passed it must be a matrix of cells
%                          with a number of rows equal to the number of nodes 
%
%
% outputs:
%
% "strFileName.NodesPositions", containing:
% * 1st column = x-values
% * 2nd column = y-values
% * from 3rd column = vector of the status
%
% "strFileName.LinksList", containing a list of rows, each row indicating a
% link as follows:
% * x_start, y_start, x_start - x_end, y_start - y_end
%
%
function ExportGraph(	strFileName,				...
						tGraph,						...
						aacNodesStatus,				... optional
						bUseSemicolonsAsSeparators	) % optional
	%
	% compute if the status of the nodes has been provided or not
	bNodesStatusIsPresent		= ( nargin > 2 );
	bInfoOnSeparatorsIsPresent	= ( nargin > 3 );
	%
	if( bNodesStatusIsPresent )
		%
		% double-check that the statuses are passed correctly
		if( size( aacNodesStatus, 1 ) ~= tGraph.iNumberOfNodes )
			%
			warning('Somebody has messed up with the dimensions of the aacNodesStatus matrix passed to MatlabToTikZ.ExportGraph(), since it has a number of rows different from the number of nodes of the network');
			iStatusDimensionality = -1;
			%
		else%
			%
			iStatusDimensionality = size( aacNodesStatus, 2 );
			%
		end;%
		%
	end;% if( bNodesStatusIsPresent )
	%
	if( ~bInfoOnSeparatorsIsPresent )
		%
		bUseSemicolonsAsSeparators = false;
		%
	end;%
	%
	%
	% create the filenames
	strNodesPositionsFileName	= strcat(strFileName, '.NodesPositions');
	strLinksListFileName		= strcat(strFileName, '.LinksList');
	%
	% ---------------------------------------------------------------------
	% open the file
	fidNodes = fopen(strNodesPositionsFileName, 'w');
	%
	% header -- notice that it depends on the existence or not of the
	% status vector
	if( bUseSemicolonsAsSeparators )
		%
		fprintf(fidNodes, 'x;y');
		%
	else%
		%
		fprintf(fidNodes, 'x y');
		%
	end;%
	%
	if( bNodesStatusIsPresent )
		%
		for iDimension = 1:iStatusDimensionality;
			%
			if( bUseSemicolonsAsSeparators )
				%
				fprintf(fidNodes, ';s%d', iDimension);
				%
			else%
				%
				fprintf(fidNodes, ' s%d', iDimension);
				%
			end;%
			%
		end;%
		%
	end;%
	%
	fprintf(fidNodes, '\n');
	%
	% write the data
	for iNode = 1:tGraph.iNumberOfNodes;
		%
		if( bUseSemicolonsAsSeparators )
			%
			fprintf(fidNodes, '%.5f;%.5f', tGraph.aafNodesPositions(iNode, 1), tGraph.aafNodesPositions(iNode, 2) );
			%
		else%
			%
			fprintf(fidNodes, '%.5f %.5f', tGraph.aafNodesPositions(iNode, 1), tGraph.aafNodesPositions(iNode, 2) );
			%
		end;%
		%
		if( bNodesStatusIsPresent )
			%
			for iDimension = 1:iStatusDimensionality;
				%
				uCellContent = aacNodesStatus{iNode, iDimension};
				if( iscell(uCellContent) )
					%
					uCellContent = uCellContent{1};
					%
				end;%
				%
				if( bUseSemicolonsAsSeparators )
					%
					fprintf(fidNodes, ';');
					%
				else%
					%
					fprintf(fidNodes, ' ');
					%
				end;%
				%
				if( ischar(uCellContent) )
					%
					fprintf(fidNodes, '%s', uCellContent);
					%
				elseif( isstring(uCellContent) )
					%
					fprintf(fidNodes, '%s', uCellContent);
					%
				elseif( isfloat(uCellContent) )
					%
					fprintf(fidNodes, '%.5f', uCellContent);
					%
				elseif( isinteger(uCellContent) )
					%
					fprintf(fidNodes, '%.5f', uCellContent);
					%
				else%
					%
					fprintf(fidNodes, '?');
					warning('Problem in the aacNodesStatus matrix passed to MatlabToTikZ.ExportGraph(): it has some cells that I don''t know how to convert in string.');
					%
				end;%
				%
			end;%
			%
		end;%
		%
		fprintf(fidNodes, '\n');
		%
	end;%
	%
	% close the file
	fclose(fidNodes);
	clear fidNodes;
	%
	fprintf('%s exported\n', strNodesPositionsFileName);
	%
	% ---------------------------------------------------------------------
	% open the file
	fidLinks = fopen(strLinksListFileName, 'w');
	%
	% header
	fprintf(fidLinks, 'xstart ystart towardsx towardsy\n');
	%
	for iNodeA = 1:tGraph.iNumberOfNodes;
		%
		for iNodeB = 1:tGraph.iNumberOfNodes;
			%
			if(iNodeA == iNodeB)
				%
				continue;
				%
			else%
				%
				% check if there is the link
				if( tGraph.aaiAdjacencyMatrix(iNodeA, iNodeB) )
					%
					fprintf(fidLinks, '%.5f %.5f %.5f %.5f',											...
							tGraph.aafNodesPositions(iNodeA, 1),										...
							tGraph.aafNodesPositions(iNodeA, 2),										...
							tGraph.aafNodesPositions(iNodeB, 1) - tGraph.aafNodesPositions(iNodeA, 1),	...
							tGraph.aafNodesPositions(iNodeB, 2) - tGraph.aafNodesPositions(iNodeA, 2));
					%
					fprintf(fidLinks, '\n');
					%
				end;%
				%
			end;% nodes are the same
			%
		end;% cycle on node B
		%
	end;% cycle on node A
	%
	% close the file
	fclose(fidLinks);
	clear fidLinks;
	%
	fprintf('%s exported\n', strLinksListFileName);
	%
end %

