LDFLAGS        = -lm
CFLAGS        = -fPIC -O2
CC             = gcc

all: RNAdegree

RNAdegree: convert_Vienna.o IL.o HP.o energy_par.o misc.o params.o utils.o ms2_Turner.o ms1_Turner.o fold_vars.o fold.o expNumNborTurner.o read_epars.o
	$(CC) convert_Vienna.o IL.o HP.o energy_par.o misc.o params.o utils.o ms2_Turner.o ms1_Turner.o fold_vars.o fold.o expNumNborTurner.o read_epars.o $(LDFLAGS) -o RNAdegree 

ms2_Turner.o: ms2_Turner.c ms2_Turner.h
	$(CC) $(CFLAGS) -c ms2_Turner.c

ms1_Turner.o: ms1_Turner.c ms1_Turner.h
	$(CC) $(CFLAGS) -c ms1_Turner.c

misc.o: misc.c misc.h
	$(CC) $(CFLAGS) -c misc.c

convert_Vienna.o: convert_Vienna.c convert_Vienna.h fold.h energy_const.h  \
	fold_vars.h pair_mat.h params.h 
	$(CC) $(CFLAGS) -c convert_Vienna.c

HP.o: HP.c HP.h convert_Vienna.h pair_mat.h
	$(CC) $(CFLAGS) -c HP.c

IL.o: IL.c IL.h
	$(CC) $(CFLAGS) -c IL.c

fold_vars.o: fold_vars.c fold_vars.h
	$(CC) $(CFLAGS) -c fold_vars.c

fold.o: fold.c utils.h energy_par.h fold_vars.h pair_mat.h params.h
	$(CC) $(CFLAGS) -c fold.c

utils.o: utils.c config.h
	$(CC) $(CFLAGS) -c utils.c

energy_par.o: energy_par.c energy_const.h 
	$(CC) $(CFLAGS) -c energy_par.c

params.o: params.c config.h energy_par.h fold_vars.h utils.h params.h
	$(CC) $(CFLAGS) -c params.c

read_epars.o: read_epars.c 
	$(CC) $(CFLAGS) -c read_epars.c
	
expNumNborTurner.o: expNumNborTurner.c 
	$(CC) $(CFLAGS) -c expNumNborTurner.c

clean:
	rm -fr *.o RNAdegree


