
/* A Bison parser, made by GNU Bison 2.4.1.  */

/* Skeleton interface for Bison's Yacc-like parsers in C
   
      Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.
   
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.
   
   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */


/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     GROUP_NAME = 258,
     STRING = 259,
     IDENTIFIER = 260,
     VARIABLE = 261,
     SCALAR = 262,
     INTEGER = 263,
     LOGICAL = 264,
     LE = 265,
     GE = 266,
     NE = 267,
     EQ = 268,
     AND = 269,
     OR = 270,
     ACOS = 271,
     ACOSH = 272,
     ASIN = 273,
     ASINH = 274,
     ATAN = 275,
     ATANH = 276,
     CBRT = 277,
     CEIL = 278,
     COS = 279,
     COSH = 280,
     ERFC = 281,
     ERF = 282,
     EXP = 283,
     EXPM1 = 284,
     FABS = 285,
     FLOOR = 286,
     J0 = 287,
     J1 = 288,
     LGAMMA = 289,
     LOG10 = 290,
     LOG1P = 291,
     LOGB = 292,
     LOG = 293,
     SIN = 294,
     SINH = 295,
     SQRT = 296,
     TAN = 297,
     TANH = 298,
     Y0 = 299,
     Y1 = 300,
     RINT = 301
   };
#endif



#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
{

/* Line 1676 of yacc.c  */
#line 388 "parse.y"
 
  int logical_type;  
  int integer_type; 
  double scalar_type;  
  char * string_type; 
  char * subgroup_type;
  struct node_expr * node_type;
  


/* Line 1676 of yacc.c  */
#line 109 "parse.tab.h"
} YYSTYPE;
# define YYSTYPE_IS_TRIVIAL 1
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
#endif

extern YYSTYPE yylval;

