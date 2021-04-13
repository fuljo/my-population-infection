# Build this image:  docker build -t my-population-infection .

FROM ubuntu:20.04 AS openmpi

ENV USER mpirun

ENV DEBIAN_FRONTEND=noninteractive \
    HOME=/home/${USER} 


RUN apt-get update -y && \
    apt-get install -y --no-install-recommends sudo apt-utils && \
    apt-get install -y --no-install-recommends openssh-server dnsutils \
        build-essential gcc gfortran binutils \
        libopenmpi-dev && \
    apt-get clean && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /var/run/sshd
RUN echo 'root:${USER}' | chpasswd
RUN sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

ENV NOTVISIBLE "in users profile"
RUN echo "export VISIBLE=now" >> /etc/profile

# ------------------------------------------------------------
# Add an 'mpirun' user
# ------------------------------------------------------------

RUN adduser --disabled-password --gecos "" ${USER} && \
    echo "${USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers

# ------------------------------------------------------------
# Set-Up SSH with previously generated keys
# ------------------------------------------------------------

ENV SSHDIR ${HOME}/.ssh/

RUN mkdir -p ${SSHDIR}

RUN echo "StrictHostKeyChecking no\n" >> ${SSHDIR}/config
ADD ssh/id_rsa ${SSHDIR}/id_rsa
ADD ssh/id_rsa.pub ${SSHDIR}/id_rsa.pub
ADD ssh/id_rsa.pub ${SSHDIR}/authorized_keys

RUN chmod -R 600 ${SSHDIR}* && \
    chown -R ${USER}:${USER} ${SSHDIR}

EXPOSE 22

# ------------------------------------------------------------
# Configure OpenMPI
# ------------------------------------------------------------

USER root

RUN rm -fr ${HOME}/.openmpi && mkdir -p ${HOME}/.openmpi
RUN chown -R ${USER}:${USER} ${HOME}/.openmpi

FROM openmpi AS my-population-infection

# ------------------------------------------------------------
# Copy and compile application code
# ------------------------------------------------------------

USER mpirun
ENV APP  my-population-infection

ENV ROOT_HOSTNAME mpi-root
ENV WORKER_HOSTNAME mpi-worker

WORKDIR /home/${USER}/

# Copy source code
COPY src/ src/

# Compile source code
RUN make -C src

# Move executable to HOME
USER root
RUN chown -R ${USER}:${USER} src/
USER mpirun
RUN mv src/${APP} ./

# ------------------------------------------------------------
# Entrypoint script
# ------------------------------------------------------------

COPY docker_entrypoint.sh ./

ENTRYPOINT ["./docker_entrypoint.sh"]
