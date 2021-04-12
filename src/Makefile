CC = mpicc
CFLAGS = -std=gnu11 -g -Wall
LDLIBS = -lm

exec = my-population-infection
objects = my-population-infection.o config.o csv.o individual.o mpi-datatypes.o world.o log.o

$(exec): $(objects)
	$(CC) $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -o $@

config.o: config.c config.h utils.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

csv.o: csv.c csv.h individual.h utils.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

individual.o: individual.c individual.h utils.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

log.o: log.c log.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -DLOG_USE_COLOR -c $< -o $@

mpi-datatypes.o: mpi-datatypes.c mpi-datatypes.h config.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

my-population-infection.o: my-population-infection.c config.h csv.h individual.h mpi-datatypes.h utils.h world.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

world.o: world.c world.h utils.h
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@

run: $(exec)
	@mpirun $(exec) \
		-N 100 \
		-I 10 \
		-W 10000 \
		-L 3000 \
		-w 5000 \
		-l 3000 \
		-v 7 \
		-d 5 \
		--t-infection=$$((10 * 60)) \
		--t-recovery=$$((60 * 60)) \
		--t-immunity=$$((3 * 60 * 60)) \
		--sim-step=1 \
		--sim-length=2

clean:
	rm -f $(exec) *.o *.gch