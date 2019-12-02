srcdir      = .
RM          = rm
CXX         = g++
CC          = g++
INCLUDES    = -Incl -I../libxml2/include -I../main
CXXFLAGS    = -O3 -ffast-math -Wall $(INCLUDES)
LDFLAGS     = -lreadline -lm 
LIBS        = -lreadline -lm 
INSTALL     = @INSTALL@
BISON       = bison
FLEX        = flex
SRCS        = grammar.y lex.l 
OBJECTS     = grammar.tab.o lex.yy.o
BISONFLAGS  = -y -v

CFLAGS = $(CXXFLAGS)

all: grammar.tab.c lex.yy.cpp rename apply clean

grammar.tab.c grammar.tab.h: $(srcdir)/grammar.y
	$(BISON) -o $(srcdir)/grammar.tab.c --defines $(BISONFLAGS) $(srcdir)/grammar.y   

lex.yy.cpp: $(srcdir)/lex.l
	$(FLEX) -o $(srcdir)/lex.yy.cpp $(srcdir)/lex.l

rename: 
	mv $(srcdir)/grammar.tab.c $(srcdir)/grammar.tab.cpp                              

apply: 
	mv $(srcdir)/grammar.tab.h   $(srcdir)/../revlanguage/parser/grammar.tab.h
	mv $(srcdir)/grammar.tab.cpp $(srcdir)/../revlanguage/parser/grammar.tab.cpp
	mv $(srcdir)/lex.yy.cpp      $(srcdir)/../revlanguage/parser/lex.yy.cpp
	rm $(srcdir)/grammar.output

clean: 
	rm -f $(srcdir)/grammar.output $(srcdir)/grammar.tab.c $(srcdir)/grammar.tab.cpp $(srcdir)/grammar.tab.h $(srcdir)/lex.yy.cpp

