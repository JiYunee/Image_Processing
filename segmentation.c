#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>

#define BYTE    unsigned char
#define maxlevel 255

#define H 243
#define W 303

void otsu() {

    FILE* infile;


    if ((infile = fopen("hand.raw", "rb")) == NULL) {
        printf("No Image File\n");
        return;
    }

    BYTE* inImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
    fread(inImg, sizeof(BYTE), W * H, infile);

    BYTE* outImg = (BYTE*)malloc(sizeof(BYTE) * W * H);

	int hist[256] = { 0 };
	float mean = 0.0;
	float var[256] = { 0.0 };

	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			hist[inImg[i * W + j]]++;

	for (int i = 0; i < 256; i++)
		mean += (float)(i * hist[i]);

	mean = mean / (float)(H * W);

	float p1 = 0, m1 = 0, m2 = 0, p = 0, m = 0;

	for (int t = 0; t < 256; t++) {
		p1 = p + (float)hist[t] / (float)(H * W);
		if (p1 == 0.) 
			m1 = 0;
		else 
			m1 = (p * m + (float)(t * hist[t]) / (float)(H * W)) / p1;
		if (p1 == 1.) 
			m2 = 0;
		else 
			m2 = (mean - p1 * m1) / (1. - p1);

		var[t] = p1 * (m1 - mean) * (m1 - mean) + (1. - p1) * (m2 - mean) * (m2 - mean);
		p = p1;
		m = m1;
	}

	int threshold = 0;
	float max = 0;

	for (int t = 0; t < 256; t++) {
		if (max < var[t]) {
			max = var[t];
			threshold = t;
		}
	}
	
	printf("==== Otsu =====\n");
	printf("¹®ÅÎ°ª : %d\n", threshold);
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++) {
			if (inImg[i * W + j] < threshold)
				outImg[i * W + j] = 0;
			else outImg[i * W + j] = 255;
		}


	FILE* out = fopen("hand_otsu.pgm", "wb");
	fprintf(out, "P5\n");
	fprintf(out, "%d %d\n", W, H);
	fprintf(out, "%d\n", maxlevel);
	fwrite(outImg, sizeof(char), H * W, out);
	fclose(out);

    fclose(infile);
    free(inImg);
    free(outImg);
}

void efficient2Pass() 
{
	FILE* infile;

	if ((infile = fopen("hand_otsu.pgm", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	int tmp = 0;
	char chtmp[5];
	fscanf(infile, "%s", chtmp);
	fscanf(infile, "%d %d", &tmp, &tmp);
	fscanf(infile, "%d", &tmp);


	BYTE* inImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* outImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* labelImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	memset(labelImg, 0, sizeof(BYTE) * W * H);

	fread(inImg, sizeof(BYTE), W * H, infile);

	// numbering
	int number = 1;
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			if (inImg[i * W + j] != 0)
				labelImg[i * W + j] = number++;

	// 2-pass  labeling
	int check = -1;
	while (check != 0)
	{
		// top-down
		for (int i = 0; i < H - 1; i++)
			for (int j = 1; j < W - 1; j++) {
				if ((labelImg[i * W + j] != 0) && (labelImg[i * W + (j + 1)] != 0) && (labelImg[i * W + j] != labelImg[i * W + (j + 1)]))
				{
					if (labelImg[i * W + j] < labelImg[i * W + (j + 1)]) labelImg[i * W + (j + 1)] = labelImg[i * W + j];
					else labelImg[i * W + j] = labelImg[i * W + (j + 1)];
				}
				if ((labelImg[i * W + j] != 0) && (labelImg[(i + 1) * W + (j + 1)] != 0) && (labelImg[i * W + j] != labelImg[(i + 1) * W + (j + 1)]))
				{
					if (labelImg[i * W + j] < labelImg[(i + 1) * W + (j + 1)]) labelImg[(i + 1) * W + (j + 1)] = labelImg[i * W + j];
					else labelImg[i * W + j] = labelImg[(i + 1) * W + (j + 1)];
				}
				if ((labelImg[i * W + j] != 0) && (labelImg[(i + 1) * W + (j - 1)] != 0) && (labelImg[i * W + j] != labelImg[(i + 1) * W + (j - 1)]))
				{
					if (labelImg[i * W + j] < labelImg[(i + 1) * W + j]) labelImg[(i + 1) * W + j] = labelImg[i * W + j];
					else labelImg[i * W + j] = labelImg[(i + 1) * W + j];
				}
				if ((labelImg[i * W + j] != 0) && (labelImg[(i + 1) * W + (j - 1)] != 0) && (labelImg[i * W + j] != labelImg[(i + 1) * W + (j - 1)]))
				{
					if (labelImg[i * W + j] < labelImg[(i + 1) * W + (j - 1)]) labelImg[(i + 1) * W + (j - 1)] = labelImg[i * W + j];
					else labelImg[i * W + j] = labelImg[(i + 1) * W + (j - 1)];
				}
			}

		//bottom-up
		for (int i = H - 1; i >= 1; i--)
			for (int j = W - 2; j >= 1; j--) {
				if ((labelImg[i * W + j] != 0) && (labelImg[i * W + (j - 1)] != 0) && (labelImg[i * W + j] != labelImg[i * W + (j - 1)]))
				{
					if (labelImg[i * W + j] < labelImg[i * W + (j - 1)]) labelImg[i * W + (j - 1)] = labelImg[i * W + j];
					else labelImg[i * W + j] = labelImg[i * W + (j - 1)];
				}
				if ((labelImg[i * W + j] != 0) && (labelImg[(i - 1) * W + (j - 1)] != 0) && (labelImg[i * W + j] != labelImg[(i - 1) * W + (j - 1)]))
				{
					if (labelImg[i * W + j] < labelImg[(i - 1) * W + (j - 1)]) labelImg[(i - 1) * W + (j - 1)] = labelImg[i * W + j];
					else labelImg[i * W + j] = labelImg[(i - 1) * W + (j - 1)];
				}
				if ((labelImg[i * W + j] != 0) && (labelImg[(i - 1) * W + j] != 0) && (labelImg[i * W + j] != labelImg[(i - 1) * W + j]))
				{
					if (labelImg[i * W + j] < labelImg[(i - 1) * W + j]) labelImg[(i - 1) * W + j] = labelImg[i * W + j];
					else labelImg[i * W + j] = labelImg[(i - 1) * W + j];
				}
				if ((labelImg[i * W + j] != 0) && (labelImg[(i - 1) * W + (j + 1)] != 0) && (labelImg[i * W + j] != labelImg[(i - 1) * W + (j + 1)]))
				{
					if (labelImg[i * W + j] < labelImg[(i - 1) * W + (j + 1)]) labelImg[(i - 1) * W + (j + 1)] = labelImg[i * W + j];
					else labelImg[i * W + j] = labelImg[(i - 1) * W + (j + 1)];
				}
			}

		check = 0;
		for (int i = 0; i < H - 1; i++)
			for (int j = 1; j < W - 1; j++) {
				if ((labelImg[i * W + j] != 0) && (labelImg[i * W + (j + 1)] != 0) && (labelImg[i * W + j] != labelImg[i * W + (j + 1)])) check++;
				if ((labelImg[i * W + j] != 0) && (labelImg[(i + 1) * W + j] != 0) && (labelImg[i * W + j] != labelImg[(i + 1) * W + j])) check++;
				if ((labelImg[i * W + j] != 0) && (labelImg[(i + 1) * W + (j + 1)] != 0) && (labelImg[i * W + j] != labelImg[(i + 1) * W + (j + 1)])) check++;
				if ((labelImg[i * W + j] != 0) && (labelImg[(i + 1) * W + (j - 1)] != 0) && (labelImg[i * W + j] != labelImg[(i + 1) * W + (j - 1)])) check++;
			}
	}

	// histogram 
	int hist[256] = { 0 };

	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			hist[labelImg[i * W + j]]++;

	int max_val = 0, max_idx = 0;

	for (int i = 1; i < 256; i++)
		if (hist[i] > max_val)
		{
			max_val = hist[i];
			max_idx = i;
		}

	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
		{
			if (labelImg[i * W + j] != max_idx)
				outImg[i * W + j] = 0;
			else
				outImg[i * W + j] = 255;
		}


	FILE* out = fopen("hand_efficient2Pass.pgm", "wb");
	fprintf(out, "P5\n");
	fprintf(out, "%d %d\n", W, H);
	fprintf(out, "%d\n", maxlevel);
	fwrite(outImg, sizeof(char), H * W, out);
	fclose(out);

	fclose(infile);
	free(inImg);
	free(outImg);
	free(labelImg);
}

void main() 
{
	otsu();
	efficient2Pass();
}
