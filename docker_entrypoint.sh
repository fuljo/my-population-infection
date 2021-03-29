#!/bin/bash
if [ -z $IS_HEAD ];
then
    # WORKER
    # Start ssh daemon
    sudo /usr/sbin/sshd -D;
else
    # HEAD
    # Start ssh daemon
    sudo /usr/sbin/sshd;
    # Produce the list of host ip addresses
    dig +short $HEAD_HOSTNAME $NODE_HOSTNAME > mpihosts;
    # Run the application (command-line arguments will be passed by docker)
    mpirun --hostfile mpihosts $MPIRUN_OPTIONS $APP $@;
fi