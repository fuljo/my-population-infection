#!/bin/bash
if [ -z $IS_ROOT ];
then
    # WORKER
    # Start ssh daemon
    sudo /usr/sbin/sshd -D;
else
    # ROOT
    # Start ssh daemon
    sudo /usr/sbin/sshd;
    # Produce the list of host ip addresses
    dig +short $ROOT_HOSTNAME $WORKER_HOSTNAME > mpihosts;
    # Run the application (command-line arguments will be passed by docker)
    mpirun --hostfile mpihosts $MPIRUN_OPTIONS $APP $@;
fi