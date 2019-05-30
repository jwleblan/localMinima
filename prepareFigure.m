% This function prepares a figure for print
function prepareFigure(hFig)

% Default Values
if ~nargin
    hFig = get(0,'currentFigure');
end
fontName = getFont('Transitional');

for hAxes = findall(hFig,'Type','axes')'
    
    % Set font of all text objects
    fontSize = 20;
    fontWeight = 'bold';
    
    if ~isempty(fontName)
        set(hAxes,'FontName',fontName);
    end
    set(hAxes,'FontSize',fontSize);
    set(hAxes,'FontWeight',fontWeight);
    
    for hText = findall(hFig,'Type','text')'
        if ~isempty(fontName)
            set(hText,'FontName',fontName);
        end
        set(hText,'FontSize',fontSize);
        set(hText,'FontWeight',fontWeight);
    end
    
    % Setup any colorbars
    for hCBar = findall(hFig,'Type','colorbar')'
        if ~isempty(fontName)
            set(hCBar,'FontName',fontName);
        end
        set(hCBar,'FontSize',fontSize);
        set(hCBar,'FontWeight',fontWeight);
    end
    
    % Set all line widths
    lineWidth = 2;
    for hLine = findobj(hAxes,'Type','Line')'
        set(hLine,'LineWidth',lineWidth);
    end
    
end

end
