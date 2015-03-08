OPTS = -module liveViz
CHARMC=/home/ralfgunter/charmworks/ref_charm/bin/charmc $(OPTS)
TESTOPTS = ++local

OBJS = ring.o

all: ring

ring: $(OBJS)
	$(CHARMC) -g -language charm++ -o ring $(OBJS)

ring.decl.h: ring.ci
	$(CHARMC)  ring.ci

clean:
	rm -f *.decl.h *.def.h conv-host *.o ring charmrun

cleanp:
	rm -f *.sts *.gz *.projrc *.topo *.out

ring.o: ring.C ring.h ring.decl.h
	$(CHARMC) -c ring.C

test: all
	./charmrun ring $(TESTOPTS) +p4
