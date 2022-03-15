all: callm_db_build callm_db_val iso_gt_mtar
	@echo "CallM build completed."

callm_db_build:  ./src/callm_db_build.cpp Makefile
	g++ -std=c++11 ./src/callm_db_build.cpp -o ./bin/callm_db_build -O3 -lpthread

callm_db_val: ./src/callm_db_val.cpp Makefile
	g++ -std=c++11 ./src/callm_db_val.cpp -o ./bin/callm_db_val -O3 -lpthread

iso_gt_mtar: ./src/callm_db_val.cpp Makefile
	g++ -std=c++11 ./src/iso_gt_mtar.cpp -o ./bin/iso_gt_mtar -O3 -lpthread

clean:
	rm ./bin/callm_db_build ./bin/callm_db_val ./bin/iso_gt_mtar
