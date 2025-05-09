\ProvidesClass{lkx_paper}[Simple Paper Class]
\NeedsTeXFormat{LaTeX2e}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions\relax
\LoadClass[12pt,letterpaper]{article}
\RequirePackage{expl3,xparse,xstring}
\ExplSyntaxOn

% Font / Layout
\usepackage[margin=1in, includefoot]{geometry}
\usepackage[T1]{fontenc}
\usepackage{mlmodern}
\usepackage{microtype}
\usepackage[indent]{parskip}
%
% Page Numbers
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyfoot[R]{\thepage}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}

% Title
\AtBeginDocument{\maketitle\thispagestyle{fancy}}

% Lists
\usepackage[shortlabels]{enumitem}
\setlist[enumerate]{topsep=1ex,itemsep=1ex,partopsep=1ex,parsep=1ex}
\setlist[itemize]{topsep=1ex,itemsep=1ex,partopsep=1ex,parsep=1ex}

% Consistent bold math
\NewCommandCopy{\oldbfseries}{\bfseries}
\renewcommand{\bfseries}{\oldbfseries\boldmath}

% Bibliography
\usepackage[style=numeric, natbib=true, isbn=false, url=false, doi=false, eprint=false]{biblatex}
\AtEveryBibitem{\clearfield{note} \clearfield{month}}
\DeclareFieldFormat{labelnumber}{\normalfont #1}
\renewcommand*\mkbibbrackets[1]{\normalfont [#1]}

\NewDocumentCommand{
	\lkxbib}{m}{\bibliography{#1}{}
	\AtEndDocument{\renewcommand\bibname{References} 
	\vspace{1em}
	\printbibliography[title=References, heading = subbibliography]}}

% General math
\usepackage{amsmath, amssymb}

% cleveref (must be loaded last)
\usepackage{bookmark}
\hypersetup{colorlinks=true,
	linkcolor=lkxemph,
	citecolor=lkxemph,
	filecolor=lkxemph,
	urlcolor=lkxemph}

\usepackage[capitalize,nameinlink,noabbrev]{cleveref}
\renewcommand{\chapterautorefname}{Chapter}
\renewcommand{\sectionautorefname}{Section}
\renewcommand{\subsectionautorefname}{Section}
\crefname{equation}{}{}
\crefformat{equation}{\normalfont{#2(#1)#3}}
