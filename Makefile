standalone:
	$(MAKE) -C src

docker: ssh/id_rsa ssh/id_rsa.pub
	docker build -t my-population-infection .

compose:
	docker-compose up \
		--abort-on-container-exit \
		--exit-code-from mpi_root
	docker-compose down

ssh/id_rsa ssh/id_rsa.pub:
	ssh-keygen -C mpiuser@mpi-root -t rsa -b 4096 -f ssh/id_rsa