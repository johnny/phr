COMPILE   = g++
LINK      = g++

CFLAGS    = -g -O0 -c
#CFLAGS    = -O3 -c
LLDFLAGS  = 
LIBS      = -lm -lpthread 

# set application name
EXAMPLE_APPL     = example
EXAMPLE_OBJECTS  = example.o basicthread.o semaphore.o condition.o mutex.o barrier.o flag.o

PETERSON_APPL     = peterson
PETERSON_OBJECTS  = peterson.o basicthread.o semaphore.o condition.o mutex.o barrier.o flag.o

CHECKPETERSON_APPL     = checkpeterson
CHECKPETERSON_OBJECTS  = checkpeterson.o basicthread.o semaphore.o condition.o mutex.o barrier.o flag.o

SCALARPRODUCT_APPL     = scalarproduct
SCALARPRODUCT_OBJECTS  = scalarproduct.o basicthread.o semaphore.o condition.o mutex.o barrier.o flag.o

RING_APPL     = ring
RING_OBJECTS  = ring.o basicthread.o semaphore.o condition.o mutex.o barrier.o flag.o

PRODCONS_APPL	 = producerconsumer
PRODCONS_OBJECTS = basicthread.o semaphore.o condition.o mutex.o barrier.o flag.o producer.o consumer.o prodcons.o

CONSUMER_APPL    = consumer
CONSUMER_OBJECTS = basicthread.o semaphore.o condition.o mutex.o barrier.o flag.o Consumer.o

all : $(EXAMPLE_APPL) $(PETERSON_APPL) $(CHECKPETERSON_APPL) $(SCALARPRODUCT_APPL) $(RING_APPL) $(PRODCONS_APPL)

$(EXAMPLE_APPL) : $(EXAMPLE_OBJECTS) Makefile
	$(LINK) $(LLDFLAGS) -o $(EXAMPLE_APPL) $(EXAMPLE_OBJECTS) $(LIBS)

$(PETERSON_APPL) : $(PETERSON_OBJECTS) Makefile
	$(LINK) $(LLDFLAGS) -o $(PETERSON_APPL) $(PETERSON_OBJECTS) $(LIBS)

$(CHECKPETERSON_APPL) : $(CHECKPETERSON_OBJECTS) Makefile
	$(LINK) $(LLDFLAGS) -o $(CHECKPETERSON_APPL) $(CHECKPETERSON_OBJECTS) $(LIBS)

$(SCALARPRODUCT_APPL) : $(SCALARPRODUCT_OBJECTS) Makefile
	$(LINK) $(LLDFLAGS) -o $(SCALARPRODUCT_APPL) $(SCALARPRODUCT_OBJECTS) $(LIBS)

$(RING_APPL) : $(RING_OBJECTS) Makefile
	$(LINK) $(LLDFLAGS) -o $(RING_APPL) $(RING_OBJECTS) $(LIBS)

$(PRODCONS_APPL) : $(PRODCONS_OBJECTS) Makefile
	$(LINK) $(LLDFLAGS) -o $(PRODCONS_APPL) $(PRODCONS_OBJECTS) $(LIBS) -g

$(CONSUMER_APPL) : $(CONSUMER_OBJECTS) Makefile
	$(LINK) $(LLDFLAGS) -o $(CONSUMER_APPL) $(CONSUMER_OBJECTS) $(LIBS)

.c.o:
	$(COMPILE) $(CFLAGS) $<

.cc.o:
	$(COMPILE) $(CFLAGS) $<

clean:
	rm *.o
	rm $(EXAMPLE_APPL) $(RING_APPL) $(PETERSON_APPL) $(CHECKPETERSON_APPL) $(SCALARPRODUCT_APPL)
