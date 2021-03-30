standalone:
	$(MAKE) -C src

docker: ssh/id_rsa ssh/id_rsa.pub
	docker build -t my-population-infection .

ssh/id_rsa ssh/id_rsa.pub:
	ssh-keygen -C mpiuser@mpi-root -t rsa -b 4096 -f ssh/id_rsa