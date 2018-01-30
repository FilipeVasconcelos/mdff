#!/bin/bash
MPIRUN=mpirun
exe=mdff.x

usage()
{
cat quench.config
cat <<EOF
  usage: $0 --serial OR --para=NPROC
EOF
}


init()
{
    echo "-----------------------------------"
    echo " Start from config  POSFF.start    "
    echo "-----------------------------------"
    if [ -f  POSFF.start ] ; then
	cp POSFF.start POSFF 
    else
	echo " ERROR: POSFF.start is missing. STOP ."
	exit 1;
    fi
}
init_random()
{
    echo "-----------------------------------"
    echo " Start from random config          "
    echo "-----------------------------------"
    random_struct.py -n 300 -i "60 30" -t "O Si" -a "11.0757 0.0 0.0" -b "0.0 11.0757 0.0" -c "0.0 0.0 11.0757"
    cp POSFF.randomPY POSFF.start

}


docp()
{
    ODIR=$1
    mkdir $ODIR
    for file in MSDFF control.F CONTFF POSFF OSZIFF log TRAJFF cmd RESTART; do
	[ -f $file ] && mv $file $ODIR/
    done
    [ -f $ODIR/CONTFF ] && cp $ODIR/CONTFF ./POSFF && cp $ODIR/RESTART RESTART
}

DoMD()
{
    echo "-----------------------------------"
    echo " $1 "
    echo "-----------------------------------"
    NAME=$7
    INPUT=$2
    cp $INPUT control.F
    echo ${CMD} > cmd 
    source cmd > log
    docp $NAME
}


if [ $# -eq 0 ]; then
    usage 1;
fi

while [ $# -gt 0 ]; do
    
    case "$1" in
	
	--*=*)
	    optarg=`echo "$1" | /bin/sed -e 's/[-_a-zA-Z0-9]*=//' `
	    ;;
	
	*)
	    optarg=
	    ;;

    esac
    
    case "$1" in
	
	--serial)
	    serial=yes
	    ;;
	
	--para=*)
	    serial=no
	    NPROC=$optarg
	    ;;

	--help)
	    usage 0
	    ;;
	
	*)
	    echo $1
	    usage 1
	    ;;
	
    esac
    shift
done

DoALL()
{
	source DOALLMD
}


#--------------------------------------------
if [ "$serial" == "no" ] ; then
    #CMD="$MPIRUN -n $NPROC numawrap --node --exe='$exe control.F'"
    CMD="$MPIRUN -n $NPROC $exe control.F"
else
    CMD="$exe control.F" 
fi
#--------------------------------------------
#--------------------------------------------
usage
#init_random
init
#--------------------------------------------
DoALL
#--------------------------------------------
echo " END "
#--------------------------------------------
exit 0;
