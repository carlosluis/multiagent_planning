# optidef

Optidef is a small library that provides a standard set of environments for writing optimization problems. 


## Features

The most important features are:

- Different alignment alternatives between the constraints, objective function and the other problem components.

- Implemenation of a short format where "minimize" is substituted by "min" and "subject to" by "s.t."

- It references optimization problem using three different policies: 
   * No equation is referenced
   * The problem is referenced with a single label
   * Each equation has an individual reference.

- Limitless number of constraints.

- Four types of optimization problems:
   * minimize
   * maximize
   * arg mini
   * arg maxi

- The objective function can be broken in several lines without compromising the alignment or the structure of the problem.

## Usage

Import the package by directly adding \usepackage{optidef} to your LaTeX document. Consult the documentation for different examples and syntax usage.


## Syntax
    
The syntax to define an optimization problem is given by:
 
        >\begin{mini#}|sizeFormat|[constraintFormat]
            {optimizationVariable}
            {objectiveFunction \label{objectiveReference}}
            {\label{problemReference}}  
            {problemResult}
            \addConstraint{LHS.1}{RHS.1 \label{Const.1}}{extraInfo.1}
            \addConstraint{LHS.2}{RHS.2 \label{Const.2}}{extraInfo.2}
            .
            .
            \addConstraint{LHS.N}{RHS.N \label{Const.N}}{extraInfo.N}
        \end{mini#}


where mini# takes any of the following values: 

 - mini\* for no referencing
 - mini! for referencing each equation 
 - mini for referencing with a single label the whole problem. 
    
The last two defined problem parameters, "\label{optimizationProblem}"" and "optimizationResult", are mandatory to allow line breaking between the 6 parameters.

After the definition of this parameters, the environment accepts the definition of an infinite number of constraints.

Finally note that \begin{mini#} can be substituted by \begin{maxi#}, \begin{argmini#} or \begin{argmaxi#}. 

## Contact for issue reporting or suggestions

E-mail: J.LagoGarcia(at)tudelft.nl

Github: https://github.com/jeslago/optidef

## Latest stable version: Optidef 2.6

CTAN: https://www.ctan.org/pkg/optidef

## Licensing

Copyright 2017 Jesus Lago

This work may be distributed and/or modified under the conditions of the LaTeX Project Public License, either version 1.3 of this license or (at your option) any later version.
The latest version of this license is in http://www.latex-project.org/lppl.txt and version 1.3 or later is part of all distributions of LaTeX version 2005/12/01 or later.

This work has the LPPL maintenance status 'maintained'. The Current Maintainer of this work is Jesus Lago.

This work consists of the file optidef.sty.