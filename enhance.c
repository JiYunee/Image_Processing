#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define WIDTHBYTES(bits) (((bits)+31)/32*4)
#define BYTE    unsigned char
#define maxlevel 255

void Q1_histogram() {
	FILE* infile;

	if ((infile = fopen("publicSquare.raw", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	int W = 1101;
	int H = 1080;
	BYTE* IpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);

	fread(IpImg, sizeof(BYTE), W * H, infile);

	int hist[256] = { 0 };
	double w[256] = { 0. };
	int t[256] = { 0 };
	int wn_w = 0;

	// calculate histogram
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			hist[IpImg[i * W + j]]++;

	for (int i = 1; i < 256; i++)
		hist[i] += hist[i - 1];

	for (int i = 0; i < 256; i++)
		w[i] = (double)hist[i] / (H * W);

	for (int i = 0; i < 256; i++) {
		while (wn_w * wn_w / (255. * 255.) - w[i] < 0)
			wn_w++;
		t[i] = wn_w;
	}

	// Equalization
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			IpImg[i * W + j] = t[IpImg[i * W + j]];


	// Image Output
	FILE* outfile = fopen("publicSquare_spec.pgm", "wb");
	fprintf(outfile, "P5\n"); //magic no
	fprintf(outfile, "%d %d\n", W, H);
	fprintf(outfile, "%d\n", maxlevel);
	fwrite(IpImg, sizeof(char), H * W, outfile);

	fclose(outfile);
	fclose(infile);
	free(IpImg);
}

void Q2_2D_meanfilter() {

	// file read
	FILE* infile;

	if ((infile = fopen("Snow.raw", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	int W = 3136;
	int H = 2199;

	BYTE* IpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* meanFilter_Img = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(IpImg, sizeof(BYTE), W * H, infile);

	// 21*21 2D filter
	int filter = 21;
	int half = filter / 2;

	// image making
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {

			if (i < half) {

				if (j < half) {
					int tmp = 0;
					int count = 0;
					for (int k = -i; k <= half; k++) {
						for (int r = -j; r <= half; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
							count++;
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / count + 0.5);
				}
				else if (j >= half && j < W - half) {
					int tmp = 0;
					int count = 0;
					for (int k = -i; k <= half; k++) {
						for (int r = -half; r <= half; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
							count++;
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / count + 0.5);
				}
				else {
					int tmp = 0;
					int count = 0;
					for (int k = -i; k <= half; k++) {
						for (int r = -half; r < W - j; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
							count++;
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / count + 0.5);
				}

			}
			else if (i >= half && i < H - half) {

				if (j < half) {
					int tmp = 0;
					int count = 0;
					for (int k = -half; k <= half; k++) {
						for (int r = -j; r <= half; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
							count++;
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / count + 0.5);
				}
				else if (j >= half && j < W - half) {
					int tmp = 0;
					for (int k = -half; k <= half; k++) {
						for (int r = -half; r <= half; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / (filter * filter) + 0.5);
				}
				else {
					int tmp = 0;
					int count = 0;
					for (int k = -half; k <= half; k++) {
						for (int r = -half; r < W - j; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
							count++;
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / count + 0.5);
				}

			}
			else {

				if (j < half) {
					int tmp = 0;
					int count = 0;
					for (int k = -half; k < H - i; k++) {
						for (int r = -j; r <= half; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
							count++;
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / count + 0.5);
				}
				else if (j >= half && j < W - half) {
					int tmp = 0;
					int count = 0;
					for (int k = -half; k < H - i; k++) {
						for (int r = -half; r <= half; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
							count++;
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / count + 0.5);
				}
				else {
					int tmp = 0;
					int count = 0;
					for (int k = -half; k < half - i; k++) {
						for (int r = -half; r < W - j; r++) {
							tmp += IpImg[(i + k) * W + (j + r)];
							count++;
						}
					}
					meanFilter_Img[i * W + j] = (int)((double)tmp / count + 0.5);
				}

			}

		}
	}

	// Image Output
	FILE* outfile = fopen("Snow_2Dmeanfilter.pgm", "wb");
	fprintf(outfile, "P5\n"); //magic no
	fprintf(outfile, "%d %d\n", W, H);
	fprintf(outfile, "%d\n", maxlevel);
	fwrite(meanFilter_Img, sizeof(char), W * H, outfile);

	fclose(outfile);
	fclose(infile);
	free(IpImg);
	free(meanFilter_Img);
}

void Q2_1D_meanfilter() {

	// file read
	FILE* infile;

	if ((infile = fopen("Snow.raw", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	int W = 3136;
	int H = 2199;

	BYTE* IpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* tmpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* meanFilter_Img = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(IpImg, sizeof(BYTE), W * H, infile);


	// 21*21 2D filter
	int filter = 21;
	int half = filter / 2;

	// horizontal & vertical 
	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			int tmp = 0;
			if (i < half) {
				for (int k = -i; k <= half; k++)
					tmp += IpImg[(i + k) * W + j];
				tmpImg[i * W + j] = (int)((double)tmp / (half + i + 1) + 0.5);
			}
			else if (i >= half && i < H - half) {
				for (int k = -half; k <= half; k++)
					tmp += IpImg[(i + k) * W + j];
				tmpImg[i * W + j] = (int)((double)tmp / filter + 0.5);
			}
			else {
				for (int k = -half; k < H - i; k++)
					tmp += IpImg[(i + k) * W + j];
				tmpImg[i * W + j] = (int)((double)tmp / (H - i + half) + 0.5);
			}
		}
	}

	for (int i = 0; i < H; i++) {
		for (int j = 0; j < W; j++) {
			int tmp = 0;
			if (i < half) {
				for (int k = -j; k <= half; k++)
					tmp += tmpImg[i * W + j + k];
				meanFilter_Img[i * W + j] = (int)((double)tmp / (half + j + 1) + 0.5);
			}
			else if (i >= half && i < W - half) {
				for (int k = -half; k <= half; k++)
					tmp += tmpImg[i * W + j + k];
				meanFilter_Img[i * W + j] = (int)((double)tmp / filter + 0.5);
			}
			else {
				for (int k = -half; k < W - j; k++)
					tmp += tmpImg[i * W + j + k];
				meanFilter_Img[i * W + j] = (int)((double)tmp / (H - i + half) + 0.5);
			}
		}
	}

	// Image Output
	FILE* outfile = fopen("Snow_1Dmeanfilter.pgm", "wb");
	fprintf(outfile, "P5\n"); //magic no
	fprintf(outfile, "%d %d\n", W, H);
	fprintf(outfile, "%d\n", maxlevel);
	fwrite(meanFilter_Img, sizeof(char), W * H, outfile);

	fclose(outfile);
	fclose(infile);
	free(IpImg);
	free(tmpImg);
	free(meanFilter_Img);
}

void Q2_3_integral() {

	// file read
	FILE* infile;

	if ((infile = fopen("Snow.raw", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	int W = 3136;
	int H = 2199;

	BYTE* IpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	int* tmpImg = (int*)malloc(sizeof(int) * W * H);
	BYTE* integral_Img = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(IpImg, sizeof(BYTE), W * H, infile);

	for (int i = 0; i < W; i++)
		tmpImg[i] = (int)IpImg[i];

	for (int i = 1; i < H; i++)
		for (int j = 0; j < W; j++)
			tmpImg[i * W + j] = tmpImg[(i - 1) * W + j] + (int)IpImg[i * W + j];

	for (int i = 0; i < H; i++)
		for (int j = 1; j < W; j++)
			tmpImg[i * W + j] += tmpImg[i * W + j - 1];

	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++) {
			int trans = (int)(((double)tmpImg[i * W + j] / tmpImg[(H - 1) * W + (W - 1)]) * 255 + 0.5);
			if (trans < 0 || trans > 255)
				integral_Img[i * W + j] = 0;
			else
				integral_Img[i * W + j] = trans;
		}

	// Image Output
	FILE* outfile = fopen("Snow_integral.pgm", "wb");
	fprintf(outfile, "P5\n"); //magic no
	fprintf(outfile, "%d %d\n", W, H);
	fprintf(outfile, "%d\n", maxlevel);
	fwrite(integral_Img, sizeof(char), W * H, outfile);

	fclose(outfile);
	fclose(infile);
	free(IpImg);
	free(tmpImg);
	free(integral_Img);

}

void Q2_4_mean_integral() {

	// file read
	FILE* infile;

	if ((infile = fopen("Snow.raw", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	int W = 3136;
	int H = 2199;

	BYTE* IpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	int* tmpImg = (int*)malloc(sizeof(int) * W * H);
	BYTE* OutImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(IpImg, sizeof(BYTE), W * H, infile);

	int filter = 21;
	int half = filter / 2;

	for (int i = 0; i < W; i++)
		tmpImg[i] = (int)IpImg[i];

	for (int i = 1; i < H; i++)
		for (int j = 0; j < W; j++)
			tmpImg[i * W + j] = tmpImg[(i - 1) * W + j] + (int)IpImg[i * W + j];

	for (int i = 0; i < H; i++)
		for (int j = 1; j < W; j++)
			tmpImg[i * W + j] += tmpImg[i * W + j - 1];

	for (int i = 0; i < H; i++)
	{
		for (int j = 0; j < W; j++)
		{
			if (i <= half)
			{
				if (j <= half)  // case 1 = d
				{
					int tmp_s = tmpImg[(i + half) * W + j + half];  //d
					OutImg[i * W + j] = (int)((double)tmp_s / ((i + 1 + half) * (j + 1 + half)) + 0.5);
				}
				else if (j > half && j < W - half)  // case 2 => d - c
				{
					int tmp_s = tmpImg[(i + half) * W + j + half] - tmpImg[(i + half) * W + j - half - 1];
					OutImg[i * W + j] = (int)((double)tmp_s / ((i + 1 + half) * filter) + 0.5);
				}
				else  //case 3 => d - c
				{
					int tmp_s = tmpImg[(i + half) * W + W - 1] - tmpImg[(i + half) * W + j - half - 1];
					OutImg[i * W + j] = (int)((double)tmp_s / ((i + 1 + half) * (W - j + half)) + 0.5);
				}
			}
			else if (i > half && i < H - half)
			{
				if (j <= half) //case 4 => d - b
				{
					int tmp_s = tmpImg[(i + half) * W + j + half] - tmpImg[(i - half - 1) * W + j + half];
					OutImg[i * W + j] = (int)((double)tmp_s / ((j + 1 + half) * filter) + 0.5);
				}
				else if (j > half && j < W - half)  // case 5 => a + d - b - c
				{
					int tmp_s = tmpImg[(i - half - 1) * W + j - half - 1] + tmpImg[(i + half) * W + j + half] - tmpImg[(i - half - 1) * W + j + half] - tmpImg[(i + half) * W + j - half - 1];  //c
					OutImg[i * W + j] = (int)((double)tmp_s / (filter * filter) + 0.5);
				}
				else  //case 6 => a + d - b - c
				{
					int tmp_s = tmpImg[(i - half - 1) * W + j - half - 1] + tmpImg[(i + half) * W + W - 1] - tmpImg[(i - half - 1) * W + W - 1] - tmpImg[(i + half) * W + j - half - 1]; //c
					OutImg[i * W + j] = (int)((double)tmp_s / ((W - j + half) * filter) + 0.5);
				}
			}
			else
			{
				if (j <= half)  //case 7 => d - b
				{
					int tmp_s = tmpImg[(H - 1) * W + j + half] - tmpImg[(i - half - 1) * W + j + half];
					OutImg[i * W + j] = (int)((double)tmp_s / ((j + 1 + half) * (H - i + half)) + 0.5);
				}
				else if (j > half && j < W - half) //case 8 => a + d - b - c
				{
					int tmp_s = tmpImg[(i - half - 1) * W + j - half - 1] + tmpImg[(H - 1) * W + j + half] - tmpImg[(i - half - 1) * W + j + half] - tmpImg[(H - 1) * W + j - half - 1];  //c
					OutImg[i * W + j] = (int)((double)tmp_s / ((H - i + half) * filter) + 0.5);
				}
				else  //case 9  => a + d - b - c
				{
					int tmp_s = tmpImg[(i - half - 1) * W + j - half - 1] + tmpImg[(H - 1) * W + W - 1] - tmpImg[(i - half - 1) * W + W - 1] - tmpImg[(H - 1) * W + j - half - 1];
					OutImg[i * W + j] = (int)((double)tmp_s / ((H - i + half) * (W - j + half)) + 0.5);
				}
			}
		}
	}

	// Image Output
	FILE* outfile = fopen("Snow_mean_integral.pgm", "wb");
	fprintf(outfile, "P5\n"); //magic no
	fprintf(outfile, "%d %d\n", W, H);
	fprintf(outfile, "%d\n", maxlevel);
	fwrite(OutImg, sizeof(char), W * H, outfile);

	fclose(outfile);
	fclose(infile);
	free(IpImg);
	free(tmpImg);
	free(OutImg);
}

void Q3_unsharp_masking() {

	// file read
	FILE* infile;

	if ((infile = fopen("Pentagon.raw", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	int W = 512;
	int H = 512;

	BYTE* IpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* OutImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* tmpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(IpImg, sizeof(BYTE), W * H, infile);

	int filter = 5;         // 필터 크기
	int half = filter / 2;
	float lam = 0.3;           // 가중치
	int* filter_value = malloc(filter * filter * sizeof(int));

	for (int i = 0; i < filter * filter; i++)
		filter_value[i] = 1;

	// smoothing
	for (int i = half; i < H - half; i++) {
		for (int j = half; j < W - half; j++) {
			int tmp_s = 0;
			for (int k = -half; k <= half; k++) {
				for (int l = -half; l <= half; l++) {
					tmp_s += filter_value[(k + half) * filter + (l + half)] * IpImg[(i + k) * W + (j + l)];

				}
			}
			// unsharp masking
			tmpImg[i * W + j] = tmp_s / (filter * filter);
		}
	}

	// unsharp masking
	for (int i = 2; i < H - 2; i++) {
		for (int j = 2; j < W - 2; j++) {
			if (i < half || i >= H - half || j < half || j >= W - half) {
				OutImg[i * W + j] = IpImg[i * W + j];
				continue;
			}
			int val = IpImg[i * W + j] + lam * (IpImg[i * W + j] - tmpImg[i * W + j]);
			if (val > 255) val = 255;
			if (val < 0) val = 0;
			OutImg[i * W + j] = val;
		}
	}

	// image output
	FILE* outfile = fopen("Pentagon_unsharp.pgm", "wb");
	fprintf(outfile, "P5\n"); //magic no
	fprintf(outfile, "%d %d\n", W, H);
	fprintf(outfile, "%d\n", maxlevel);
	fwrite(OutImg, 1, W * H, outfile);

	fclose(outfile);
	fclose(infile);
	free(IpImg);
	free(OutImg);
}

void Q4_interpolation() {

	// file read
	FILE* infile;

	if ((infile = fopen("Lena.raw", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	int H = 512;
	int W = 512;

	int Oct = 1024;

	BYTE* IpImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* zeroOc = (BYTE*)malloc(sizeof(BYTE) * Oct * Oct);
	BYTE* firstOc = (BYTE*)malloc(sizeof(BYTE) * Oct * Oct);

	fread(IpImg, sizeof(BYTE), W * H, infile);

	// zero
	for (int i = 0; i < Oct; i++)
		for (int j = 0; j < Oct; j++) {
			zeroOc[i * Oct + j] = IpImg[(i * H / Oct) * W + (j * W / Oct)];
		}

	// first
	for (int i = 0; i < 2 * H - 1; i += 2) {
		for (int j = 0; j < 2 * W - 1; j += 2) {
			int a = IpImg[(i / 2) * W + (j / 2)];
			int b = IpImg[(i / 2) * W + (j / 2) + 1];
			int c = IpImg[((i / 2) + 1) * W + (j / 2)];
			int d = IpImg[((i / 2) + 1) * W + (j / 2) + 1];
			firstOc[(i * 1024) + j] = a;
			firstOc[(i * 1024) + j + 1] = (BYTE)(((float)a + (float)b) / 2.0f);
			firstOc[((i + 1) * 1024) + j] = (BYTE)(((float)a + (float)c) / 2.0f);
			firstOc[((i + 1) * 1024) + j + 1] = (BYTE)(((float)a + (float)b + (float)c + (float)d) / 4.0f);
		}
	}

	// image output
	FILE* zero_out = fopen("Lena_zero.pgm", "wb");
	fprintf(zero_out, "P5\n"); //magic no
	fprintf(zero_out, "%d %d\n", Oct, Oct);
	fprintf(zero_out, "%d\n", maxlevel);
	fwrite(zeroOc, 1, sizeof(BYTE) * Oct * Oct, zero_out);

	// image output
	FILE* first_out = fopen("Lena_first.pgm", "wb");
	fprintf(first_out, "P5\n"); //magic no
	fprintf(first_out, "%d %d\n", W * 2, H * 2);
	fprintf(first_out, "%d\n", maxlevel);
	fwrite(firstOc, 1, W * 2 * H * 2, first_out);

	fclose(first_out);
	fclose(zero_out);
	fclose(infile);
	free(IpImg);
	free(zeroOc);
	free(firstOc);
}

int main()
{
	Q1_histogram();

	clock_t start1 = clock();
	Q2_2D_meanfilter();
	printf("2-1) 수행시간 : %d \n", clock() - start1);

	clock_t start2 = clock();
	Q2_1D_meanfilter();
	printf("2-2) 수행시간 : %d \n", clock() - start2);

	Q2_3_integral();

	clock_t start3 = clock();
	Q2_4_mean_integral();
	printf("2-4) 수행시간 : %d \n", clock() - start3);

	Q3_unsharp_masking();
	Q4_interpolation();
}