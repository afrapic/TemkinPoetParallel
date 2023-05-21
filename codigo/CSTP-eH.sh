#
#!/bin/bash
#$ -S /bin/bash
# ... Ask for memory
#
#$ -l vf=750M
#
# ... Change to working directory
#
#$ -cwd
#
# export environment variables
#
#$ -V
#
# merge error and output files
#
#$ -j y
#
# print date and time, optional but useful
date
# 
# Set the Parallel Environment and number of procs.
#
#
#$ -pe mpi 9
export PATH=/opt/intel/impi/3.1/bin64/:$PATH
#
#	NSLOTS is the number of processes (or processors) needed
#	       set before with -pe mpi
#	NODES  is the number of computers we need to boot
#	       Usually NSLOTS > NODES in multicore  
#
#
echo '=================================='
echo '          PROCESSORS'   
echo $TMPDIR
cat $TMPDIR/machines
cat $TMPDIR/machines >| uniq.nodes
#
#	 Start the virtual machine on each NODE with ssh transport
#
#
#main_node=`hostname`
#echo "SGE script is run on $main_node"
echo "Nodes =" $NSLOTS
#mpdboot --totalnum=$NODES --mpd=/opt/intel/impi/3.1/bin64/mpd \
#           --file=uniq.nodes --rsh=ssh 
# Optional, show the participating nodes:
# 
#	Now, run the MPI-program
#
#------------------------------------------------
#mpirun -f $TMPDIR/machines -r ssh -machinefile $TMPDIR/machines -genv I_MPI_PIN disable -n $NSLOTS ./ptp
#mpirun -f $TMPDIR/machines -r ssh -machinefile $TMPDIR/machines -genv I_MPI_DEVICE ssm \
#                      -genv I_MPI_PIN disable -n $NSLOTS ./ptp
mpirun -f uniq.nodes -r ssh -machinefile uniq.nodes -genv I_MPI_DEVICE ssm  \
                      -genv I_MPI_PIN disable -genv I_MPI_DEBUG 5 -n $NSLOTS ./cstp
#------------------------------------------------
#echo '=================================='
#echo '             NODES'
#mpdtrace -l
#echo '=================================='
# print date and time again
date
