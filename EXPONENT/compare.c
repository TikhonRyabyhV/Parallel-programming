#include <stdio.h>
#include <stdlib.h>

int main () {

	FILE* real_exp = fopen("real_exp.txt", "r");
	FILE*  res_exp = fopen( "res_exp.txt", "r");

	if(real_exp == NULL || res_exp == NULL)
		return 0;


	int item = fgetc(res_exp);

	while (item != '.'){
		item = fgetc(res_exp);
	}

	item = fgetc(res_exp);

	int sample = fgetc(real_exp);

	while (sample != '.'){
		sample = fgetc(real_exp);
	}

	sample = fgetc(real_exp);

	int count = 0, num  = 0;

	while (('1' <= item && item <= '9') || item == '0'){
		if(item == sample) {
			++count;
		}

		else {
			++num;
			break;
		}

		++num;

		item   = fgetc( res_exp);
		sample = fgetc(real_exp);
	}
	
	while (('1' <= item && item <= '9') || item == '0'){
		++num;

		item   = fgetc( res_exp);
		sample = fgetc(real_exp);
	}

	fclose(real_exp);
	fclose( res_exp);

	printf("First %d digits are correct (of all %d).\n", count, num);

	return 0;

}
