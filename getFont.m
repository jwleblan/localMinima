% getFont - Font selection helper
%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   fontName = getFont();
%   fontName = getFont(name);
%   getFont(dataBaseFile);
%
% PURPOSE:
%   This function provides a consitent interface to font selection across
%   platforms by searching for available fonts within a standard font
%   classification paradigm.  The input argument can be either a font name
%   or a font class.  If a class is provided, the highest rated avialable font
%   from that class is chosen.
%
%   The available font classes are:
%   
%   --- Serif Types ---
%   
%   Transitional:
%       A mid 18th century style made popular by John Baskerville, these
%       typefaces represent the transition between old style and neoclassical.
%       The strokes normally have verticle stresses with a large contrast
%       between thick and thin strokes.  The serifs tend to be relatively thin
%       due to a preference for horizontal attacks.
%           c.f.  Baskerville, Bookman, Cheltenham, Times New Roman
%
%   Old Style:
%       This font style dominated print prior to Baskerville's influence on
%       Benjamin Franklin and printing in the west.  The weight stress is at
%       roughly 10 and 2, and there is little variation in stroke thickness.
%       Head serifs tend to be more strongly angled leading to thicker serifs
%       overall.
%           c.f  Bembo, Caslon, Centaur, Garamond
%
%   Neoclassical:
%       This late 18th centruy style was championed by Giambattista Bodoni.  The
%       verticle stresses of transitional fonts remains, but the use of weight
%       change tends to be abrupt with serifs that are rounded and sometimes
%       ball shaped.  These accents are indicative of hand printing with ink.
%           c.f.  Bondoni, Diote, Marconi
%   
%   Slab:
%       This 19th century fonts are get there "slabiness" from heavy serifs with
%       little to no bracketing.  This, in conjunction with little to no weight
%       variation make for a very mechanistic looking font.
%           c.f.  Courier, Egyptienne, Officina, Rockwell
%
%   Glyphic:
%       These typefaces have minimal weight variation and verticle stroke
%       stresses in conjunction with a characteristic triangular serifs.  The
%       overall look tends be that of a chiseled inscription.
%           c.f.  Albertus
%
%   --- Sans-Serif Types ---
%   
%   Humanist:
%       Based on the proportions of roman inscriptions, these typefaces are
%       nearest serif fonts in due to their calligraphic influence.  Some even
%       classify them as serif fonts with very short serifs and thick brackts.
%       These fonts are known for their high readability.
%           c.f.  Frutiger, Gill Sans, Optima, Thesis, Lucidia, Meta
%
%   Grotesque:
%       A mid 19th centrury style, these fonts have the greatest variation of
%       stroke width within the sans-serif family.  The "g" often gets a bowl
%       on the end of the stroke, and/or a spur to the right.  The word gothic
%       often appears in the fonts name, indicative of the styles overall
%       appearence.
%           c.f. Helvetica, Akzidenz Grotesk, Franklin Gothic
%
%   Geometric:
%       Character shapes are based on simple geometric constructions.  The
%       repetition of similar geometric patterns creates a uniform look, if not
%       somewhat less elegant than more customized styles.
%           c.f.  Avenir, Futura, Eurostile, Harmonic Sans
%
%   Square:
%       These fonts take after the grotesque style but incorporate a geometric
%       squaring of the curves.  Unlike geometric, liberty is often taken with
%       letter spacing to mimic an overall serif appearence.  These fonts tend
%       to not do well outside of a layout engine becasue of the need to adjust
%       kerning to make things work visually.
%           c.f.  Cachet
%   
% INPUT:
%   name        - Name of a font or font-class
% 
% OUTPUT:
%   fontName    - Font name available on your system
%
%
% NOTES:
%   Other styles less well suited for technical work have been omitted.
%   c.f. Vox-ATypI Classification
%
%   For more information try:   Bringhurst's "The Elements of Typographic Style"
%
%   As with most things, Java has managed to entirely screw up their
%   antialiasing as of 2014.  As a consequence, fonts that should be beautiful
%   will render terribly in Matlab (c.f. Baskerville).  The default database
%   tries to work around this a bit through its default ordering.
%   
% ASSUMPTIONS: 
%   All input variables are of the correct type, valid(if applicable),
%   and given in the correct order. 
%
%-------------------------------------------------------------------------------
%}
function fontName = getFont(str)
persistent db
persistent allFonts
persistent allFontsLow
if isempty(db)
    db = lowerDatabase( generateDataBase() );
end

if isempty(allFonts)
    % Load system fonts
    allFonts = listfonts;
    allFontsLow = lower(allFonts);
end


% Load an alternate database if needed
if nargout==0 && nargin && exist(str,'file')
    load(str,'db');
    db = lowerDatabase(db);
    return
end
numClasses = length(db);


% Get the class and font indicies
classInd = 1;
fontInd = 1;
if nargin
    str = lower(str);
    
    % Check to see if an available fontname was given directly
    fontName = allFonts(strcmp(str,allFonts));
    if ~isempty(fontName)
        return
    end
    
    % Get the font class and name index
    classInd = find(strcmp(str,{db.class}));
    if isempty(classInd)
        % Check to see if a non-system font was given
        for i = 1:numClasses
            fontInd = find(strcmp(str,db(i).fonts));
            if ~isempty(fontInd)
                classInd = i;
                break;
            end
        end
        
        if isempty(fontInd)
            error('getFont:Unknown','Unknown Font');
        end
    end
end


% Load the closest available font
for classInd = mod(classInd-1:classInd+numClasses-2 , numClasses)+1
    for fontInd = fontInd:numel(db(classInd).fonts)
        fontName = allFonts(strcmp(db(classInd).fonts{fontInd},...
            allFontsLow));
        if ~isempty(fontName)
            fontName = fontName{1};
            return;
        end
    end
end

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   db = generateDataBase;
%
% PURPOSE:
%   This function generates the default database.
%
% INPUT:
%   None 
%
% OUTPUT:
%   db      - Database structure
%
%-------------------------------------------------------------------------------
%}
function db = generateDataBase

% Generate the database (order matters)
db(1).class = 'Transitional';
db(1).fonts = {...
    'Times New Roman',...
    'Century Schoolbook',...
    'Baskerville Old Face',...
    'Georgia',...
    'Arnhem',...
    'Joanna',...
    'Caslon',...
    'Mrs Eaves',...
    'Aurora',...
    'Baskerville',...
    'Times',...
    'Bell',...
    'Bookman',...
    'Bulmer',...
    'Caledonia',...
    'Cambria',...
    'Century',...
    'Century Gothic',...
    'Clearface',...
    'Corona',...
    'Excelsior',...
    'Imprint',...
    'Literaturnaya',...
    'Miller',...
    'Monticello',...
    'New York',...
    'Perpetua',...
    'Plantin',...
    'Utopia'};
db(2).class = 'Humanist';
db(2).fonts = {...
    'Frutiger',...
    'Palatino',...      % Right in the middle
    'Optima',...
    'Stone Sans Sem ITC TT',...
    'Skia',...
    'Tahoma',...
    'Lucidia Fax',...
    'Lucidia Bright',...
    'Lucida Sans',...
    'Gill Sans',...
    'Thesis',...
    'Meta',...
    'Trinité',...
    'Lucida',...
    'Stone',...
    'Myriad',...
    'Calibri',...
    'Candara',...
    'Centaur',...
    'Bembo',...
    'Segoe',...
    'Verdana',...
    'Clearview',...
    'Lucida Grande'};
db(3).class = 'Neoclassical';
db(3).fonts = {...
    'Didot',...
    'Walbaum',...
    'Bodoni',...
    'Marconi',...
    'Swift',...
    'Base',...
    'Matrix',...
    'Surveyor',...
    'Bell MT'};
db(4).class = 'Old Style';
db(4).fonts = {...
    'Palatino',...      % Right in the middle
    'Garamond',...
    'Bembo',...
    'Sabon',...
    'Minion',...
    'Lexicon',...
    'Cooper Black',...
    'Proforma',...
    'Janson'};
db(5).class = 'Slab';
db(5).fonts = {...
    'Rockwell',...
    'Officina',...
    'Clarendon',...
    'Courier',...
    'Courier New',...
    'Egyptienne'};
db(6).class = 'Glyphic';
db(6).fonts = {...
    'Albertus',...
    'Cartier Book'};
db(7).class = 'Grotesque';
db(7).fonts = {...
    'Akzidenz Grotesk',...
    'Franklin Gothic Medium',...
    'Franklin Gothic Book',...
    'Franklin Gothic',...
    'Inconsolata',...
    'Amplitude',...
    'Bell Gothic',...
    'Trade Gothic',...
    'News Gothic',...
    'Helvetica'};
db(8).class = 'Geometric';
db(8).fonts = {...
    'Futura',...
    'Avant Garde',...
    'Eurostile',...
    'Harmonic Sans',...
    'Gotham',...
    'Avenir',...
    'Neutraface',...
    'Nobel',...
    'Bank Gothic',...
    'Kabel',...
    'Chalet'};
db(9).class = 'Square';
db(9).fonts = {...
    'Cachet'};

end


%{ 
%-------------------------------------------------------------------------------
% SYNTAX:
%   db = lowerDatabase(db);
%
% PURPOSE:
%   This function converts all fields in the database to lower
%
% INPUT:
%   db      - Database structure
%
% OUTPUT:
%   db      - Database structure
%
%-------------------------------------------------------------------------------
%}
function db = lowerDatabase(db);

for i = 1:length(db)
    db(i).class = lower(db(i).class);
    db(i).fonts = lower(db(i).fonts);
end

end

