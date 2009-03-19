include ../Makefile.in

BINDIR = ..

INCLUDES = -I../Lib 
CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)

LIBS = -lmetis -lm
LD = $(CC) $(LDOPTIONS) -L. -L..
METISLIB = ../libmetis.a


GRAPHCHKOBJ = graphchk.o io.o 

.c.o:
	$(CC) $(CFLAGS) -c $*.c

default: $(BINDIR)/graphchk$(VERNUM)

$(BINDIR)/graphchk$(VERNUM): $(GRAPHCHKOBJ) $(METISLIB)
	$(LD)  -o $@ $(GRAPHCHKOBJ) $(LIBS)
	chmod 755 $@


clean:
	rm -f *.o

realclean:
	rm -f *.o ;\
        rm -f $(BINDIR)/graphchk$(VERNUM)
