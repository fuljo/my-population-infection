CC = mpicc
CFLAGS = -g -Wall

exec = my-population-infection
objects = my-population-infection.o config.o

$(exec): $(objects)
	$(CC) $(LDFLAGS) $^ $(LOADLIBES) $(LDLIBS) -o $@

%.o: src/%.c
	$(CC) -c $(CFLAGS) $(CPPFLAGS) $< -o $@

run: $(exec)
	mpirun $(exec) \
		-N 100 \
		-I 10 \
		-W 10000 \
		-L 3000 \
		-w 5000 \
		-l 1000 \
		-v 7 \
		-d 5 \
		--t-infection=$$((10 * 60)) \
		--t-recovery=$$((60 * 60)) \
		--t-immunity=$$((3 * 60 * 60)) \
		--sim-step=60 \
		--sim-length=2

clean:
	rm -f $(exec) *.o *.gch