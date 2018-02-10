GXX = g++
#FLAG = -03

HEADS = Koumura/Model.h Koumura/Param.h Koumura/Molecule.h Koumura/Reaction.h Koumura/Result.h Koumura/Simulation.h Koumura/Stimulus.h Koumura/UnitConversion.h General/Collection/Container2.h General/General/String.h General/Random/dSFMT.h General/Random/Random.h  General/Random/dSFMT-params.h  General/Random/dSFMT-params19937.h
OBJS = Koumura/Model.o Koumura/main.o Koumura/Param.o Koumura/Molecule.o Koumura/Reaction.o Koumura/Result.o Koumura/Simulation.o Koumura/Stimulus.o Koumura/UnitConversion.o General/General/String.o General/Random/dSFMT.o General/Random/Random.o

.SUFFIXES: .cpp .o .h

../exe: $(OBJS)
	$(GXX) $(FLAG) -o $@ $^
.cpp.o:
	$(GXX) $(FLAG) -o $@ -c $<
$(OBJS): $(HEADS)
