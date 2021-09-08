standalone:
	$(MAKE) -C src

docker: ssh/id_rsa ssh/id_rsa.pub
	docker build -t fuljo/my-population-infection .

compose:
	@mkdir -p results
	@docker compose up \
		--abort-on-container-exit \
		--exit-code-from mpi-root
	@docker compose down

ssh/id_rsa ssh/id_rsa.pub:
	mkdir -p ssh
	ssh-keygen -C mpiuser@mpi-root -t rsa -b 4096 -f ssh/id_rsa