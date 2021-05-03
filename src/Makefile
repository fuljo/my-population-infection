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
	@mpirun -np 2 $(exec) \
		-N 10000 \
		-I 100 \
		-W 2000 \
		-L 4000 \
		-w 2000 \
		-l 2000 \
		-v 1.4 \
		-d 2 \
		--t-infection=$$((1 * 60)) \
		--t-recovery=$$((1 * 3600)) \
		--t-immunity=$$((4 * 3600)) \
		--sim-step=60 \
		--sim-length=1 \
		--write-trace \
		--log-level INFO

clean:
	rm -f $(exec) *.o *.gch