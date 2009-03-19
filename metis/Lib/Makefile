SHELL:=/bin/sh

SRC =  coarsen.c fm.c initpart.c match.c ccgraph.c memory.c \
       pmetis.c pqueue.c refine.c util.c timing.c debug.c \
       bucketsort.c graph.c stat.c kmetis.c kwayrefine.c \
       kwayfm.c balance.c ometis.c srefine.c sfm.c separator.c \
       mincover.c mmd.c mesh.c meshpart.c frename.c fortran.c \
       myqsort.c compress.c parmetis.c estmem.c \
       mpmetis.c mcoarsen.c mmatch.c minitpart.c mbalance.c \
       mrefine.c mutil.c mfm.c mkmetis.c mkwayrefine.c mkwayfmh.c \
       mrefine2.c minitpart2.c mbalance2.c mfm2.c \
       kvmetis.c kwayvolrefine.c kwayvolfm.c subdomains.c

OBJS:= $(patsubst %.c, $(O_DIR)%.o, $(SRC) )

$(O_DIR)%.o  : %.c
	$(CC) -c $(CFLAGS) -o $@ $<

$(LIBMETIS): $(OBJS)
	ar $(ARFLAGS) $@ $(OBJS)
	ar -ts $@

clean:
	rm -f *.o

$(O_DIR)coarsen.o : coarsen.c 
$(O_DIR)fm.o : fm.c  
$(O_DIR)initpart.o : initpart.c
$(O_DIR)match.o : match.c 
$(O_DIR)ccgraph.o : ccgraph.c 
$(O_DIR)memory.o : memory.c 
$(O_DIR)pmetis.o : pmetis.c 
$(O_DIR)pqueue.o : pqueue.c 
$(O_DIR)refine.o : refine.c 
$(O_DIR)util.o : util.c 
$(O_DIR)timing.o : timing.c 
$(O_DIR)debug.o : debug.c 
$(O_DIR)bucketsort.o : bucketsort.c 
$(O_DIR)graph.o : graph.c 
$(O_DIR)stat.o : stat.c 
$(O_DIR)kmetis.o : kmetis.c 
$(O_DIR)kwayrefine.o : kwayrefine.c 
$(O_DIR)kwayfm.o : kwayfm.c 
$(O_DIR)balance.o : balance.c 
$(O_DIR)ometis.o : ometis.c 
$(O_DIR)srefine.o : srefine.c 
$(O_DIR)sfm.o : sfm.c 
$(O_DIR)separator.o : separator.c 
$(O_DIR)mincover.o : mincover.c 
$(O_DIR)mmd.o : mmd.c 
$(O_DIR)mesh.o : mesh.c 
$(O_DIR)meshpart.o : meshpart.c 
$(O_DIR)frename.o : frename.c 
$(O_DIR)fortran.o : fortran.c 
$(O_DIR)myqsort.o : myqsort.c
$(O_DIR)compress.o : compress.c 
$(O_DIR)parmetis.o : parmetis.c 
$(O_DIR)estmem.o : estmem.c 
$(O_DIR)mpmetis.o : mpmetis.c
$(O_DIR)mcoarsen.o : mcoarsen.c 
$(O_DIR)mmatch.o : mmatch.c 
$(O_DIR)minitpart.o : minitpart.c 
$(O_DIR)mbalance.o : mbalance.c 
$(O_DIR)mrefine.o : mrefine.c 
$(O_DIR)mutil.o : mutil.c 
$(O_DIR)mfm.o : mfm.c 
$(O_DIR)mkmetis.o : mkmetis.c 
$(O_DIR)mkwayrefine.o : mkwayrefine.c 
$(O_DIR)mkwayfmh.o : mkwayfmh.c 
$(O_DIR)mrefine2.o : mrefine2.c 
$(O_DIR)minitpart2.o : minitpart2.c 
$(O_DIR)mbalance2.o : mbalance2.c
$(O_DIR)mfm2.o : mfm2.c 
$(O_DIR)kvmetis.o : kvmetis.c 
$(O_DIR)kwayvolrefine.o : kwayvolrefine.c 
$(O_DIR)kwayvolfm.o : kwayvolfm.c 
$(O_DIR)subdomains.o : subdomains.c
