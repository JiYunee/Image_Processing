#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>


#define BYTE    unsigned char
#define maxlevel 255

//Lena
#define H 512
#define W 512

//Ctest
#define CH 550
#define CW 550

#define window_size 5


void Q1_Sobel() {
	FILE* infile;

	if ((infile = fopen("Lena.pgm", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}


	BYTE* inImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* outImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(inImg, sizeof(BYTE), W * H, infile);

	int Gx, Gy, G;

	int kernel_x[3][3] = { {-1, 0, 1}, {-2, 0, 2}, {-1, 0, 1} };
	int kernel_y[3][3] = { {-1, -2, -1}, {0, 0, 0}, {1, 2, 1} };

	for (int i = 1; i < H - 1; i++) {
		for (int j = 1; j < W - 1; j++) {
			// x 방향 sobel 연산자
			Gx = kernel_x[0][0] * inImg[(i - 1) * W + (j - 1)] + kernel_x[0][1] * inImg[(i - 1) * W + j] 
				+ kernel_x[0][2] * inImg[(i - 1) * W + (j + 1)] + kernel_x[1][0] * inImg[i * W + (j - 1)]
				+ kernel_x[1][1] * inImg[i * W + j] + kernel_x[1][2] * inImg[i * W + (j + 1)]
				+ kernel_x[2][0] * inImg[(i + 1) * W + (j - 1)] + kernel_x[2][1] * inImg[(i + 1) * W + j]
				+ kernel_x[2][2] * inImg[(i + 1) * W + (j + 1)];

			// y 방향 sobel 연산자
			Gy = kernel_y[0][0] * inImg[(i - 1) * W + (j - 1)] + kernel_y[0][1] * inImg[(i - 1) * W + j]
				+ kernel_y[0][2] * inImg[(i - 1) * W + (j + 1)] + kernel_y[1][0] * inImg[i * W + (j - 1)]
				+ kernel_y[1][1] * inImg[i * W + j] + kernel_y[1][2] * inImg[i * W + (j + 1)]
				+ kernel_y[2][0] * inImg[(i + 1) * W + (j - 1)] + kernel_y[2][1] * inImg[(i + 1) * W + j]
				+ kernel_y[2][2] * inImg[(i + 1) * W + (j + 1)];

			// gradient 계산
			G = sqrt(Gx * Gx + Gy * Gy);

			// [0, 255] 정규화
			G = G > 255 ? 255 : G;
			G = G < 0 ? 0 : G;

			outImg[i * W + j] = (BYTE)G;
		}
	}

	// 파일 출력
	FILE* out = fopen("Lena_Sobel.pgm", "wb");
	fprintf(out, "P5\n"); //magic no
	fprintf(out, "%d %d\n", H, W);
	fprintf(out, "%d\n", maxlevel);
	fwrite(outImg, 1, sizeof(BYTE) * H * W, out);

	fclose(out);
	fclose(infile);
	free(inImg);
	free(outImg);
}

void Gradient_threshold(BYTE* Img) {
	int hist[256] = { 0 };
	int half_size = window_size / 2;
	for (int i = half_size; i < H - half_size; i++)
		for (int j = half_size; j < W - half_size; j++)
			hist[Img[i * W + j]]++;

	double upval = (H - 2 * half_size) * (W - 2 * half_size) * 0.30;

	int cnt = 0;
	int thres = 255;
	while (cnt < upval)
	{
		cnt += hist[thres];
		thres--;
	}
	
	for (int i=0; i<H;i++)
		for (int j = 0; j < W; j++) {
			if (Img[i * W + j] > thres)
				Img[i * W + j] = 0;
			else
				Img[i * W + j] = 255;
		}
}

void Q1_NonlinearGradient() {
	FILE* infile;

	if ((infile = fopen("Lena.pgm", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	BYTE* inImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(inImg, sizeof(BYTE), W * H, infile);

	BYTE* outImg1 = (BYTE*)malloc(sizeof(BYTE) * W * H);
	BYTE* outImg2 = (BYTE*)malloc(sizeof(BYTE) * W * H);

	int count;
	int half_size = window_size / 2;
	int threshold = 30;

	//min
	for (int i = half_size; i < H - half_size; i++) {
		for (int j = half_size; j < W - half_size; j++) {
			BYTE minvalue = 255;
			for (int k = -half_size; k <= half_size; k++)
				for (int l = -half_size; l <= half_size; l++)
					minvalue = min(minvalue, inImg[(i + k) * W + j + l]);
			outImg1[i * W + j] = inImg[i * W + j] - minvalue;
		}
	}

	Gradient_threshold(outImg1);

	//max
	for (int i = half_size; i < H - half_size; i++) {
		for (int j = half_size; j < W - half_size; j++) {
			BYTE maxvalue = 0;
			for (int k = -half_size; k <= half_size; k++) 
				for (int l = -half_size; l <= half_size; l++) 
					maxvalue = max(maxvalue, inImg[(i + k) * W + j + l]);
			outImg2[i * W + j] = maxvalue - inImg[i * W + j];
		}
	}
	Gradient_threshold(outImg2);

	// 파일 출력 min
	FILE* out1 = fopen("Lena_NonlinearGradient_min.pgm", "wb");
	fprintf(out1, "P5\n"); //magic no
	fprintf(out1, "%d %d\n", H, W);
	fprintf(out1, "%d\n", maxlevel);
	fwrite(outImg1, 1, sizeof(BYTE) * H * W, out1);

	// 파일 출력 max
	FILE* out2 = fopen("Lena_NonlinearGradient_max.pgm", "wb");
	fprintf(out2, "P5\n"); //magic no
	fprintf(out2, "%d %d\n", H, W);
	fprintf(out2, "%d\n", maxlevel);
	fwrite(outImg2, 1, sizeof(BYTE) * H * W, out2);

	fclose(out1);
	fclose(out2);
	fclose(infile);
	free(inImg);
	free(outImg1);
	free(outImg2);
}


double integral_av(int* integral, int i, int j) {
	int sum, av_temp;
	int half_size = window_size / 2;

	int i0 = i - half_size;
	int j0 = j - half_size;
	int i1 = i + half_size;
	int j1 = j + half_size;

	if (i1 > H || j1 > W)
		return 0.0;

	if (i0 <= 0)
	{
		if (j0 <= 0) {
			sum = integral[i1 * W + j1];
			av_temp = (i1 + 1) * (j1 + 1);
		}
		else if (j0 > 0 && j1 < W) {
			sum = integral[i1 * W + j1] - integral[i1 * W + j0 - 1];
			av_temp = (i1 + 1) * (j1 - j0 + 1);
		}
		else {
			sum = integral[i1 * W + W - 1] - integral[i1 * W + j0 - 1];
			av_temp = (i1 + 1) * (W - j0 + 1);
		}
	}
	else if (i0 > 0 && i1 < H)
	{
		if (j0 <= 0) {
			sum = integral[i1 * W + j1] - integral[(i0 - 1) * W + j1];
			av_temp = (i1 - i0 + 1) * (j1 + 1);
		}
		else if (j0 > 0 && j1 < W) {
			sum = integral[(i0 - 1) * W + j0 - 1] + integral[i1 * W + j1] - integral[(i0 - 1) * W + j1] - integral[i1 * W + j0 - 1];
			av_temp = (i1 - i0 + 1) * (j1 - j0 + 1);
		}
		else {
			sum = integral[(i0 - 1) * W + j0 - 1] + integral[i1 * W + W - 1] - integral[(i0 - 1) * W + W - 1] - integral[i1 * W + j0 - 1];
			av_temp = ((i1 - i0 + 1) * (W - j0 + 1));
		}
	}
	else
	{
		if (j0 <= 0) {
			sum = integral[(H - 1) * W + j1] - integral[(i0 - 1) * W + j1];
			av_temp = ((H - i0 + 1) * (j1 + 1));
		}
		else if (j0 > 0 && j1 < W){
			sum = integral[(i0 - 1) * W + j0 - 1] + integral[(H - 1) * W + j1] - integral[(i0 - 1) * W + j1] - integral[(H - 1) * W + j0 - 1];
			av_temp = ((H - i0 + 1) * (j1 - j0 + 1));
		}
		else{
			sum = integral[(i0 - 1) * W + j0 - 1] + integral[(H - 1) * W + W - 1] - integral[(i0 - 1) * W + W - 1] - integral[(H - 1) * W + j0 - 1];
			av_temp = (H - i0 + 1) * (W - j0 + 1);
		}
	}

	return (double)sum / av_temp;
}

void Laplacian_threshold(BYTE* Img, int* tempImg, double* l_value) {
	int half_size = window_size / 2;
	double g_value = 0.0;

	int* integral = (int*)malloc(sizeof(int) * W * H);
	double* l_av = (double*)malloc(sizeof(double) * W * H);

	// integral image
	for (int i = 0; i < W; i++)
		integral[i] = tempImg[i];

	for (int i = 1; i < H; i++)
		for (int j = 0; j < W; j++)
			integral[i * W + j] = integral[(i - 1) * W + j] + tempImg[i * W + j];

	for (int i = 0; i < H; i++)
		for (int j = 1; j < W; j++)
			integral[i * W + j] += integral[i * W + j - 1];

	// integral image average
	for (int i=0;i<H;i++)
		for (int j = 0; j < W; j++) 
			l_av[i * W + j] = integral_av(integral, i, j);

	for (int i = half_size; i < H - half_size; i++)
		for (int j = half_size; j < W - half_size; j++) {
			double sum_temp = 0.0;
			for (int k = -half_size; k <= half_size; k++)
				for (int l = -half_size; l <= half_size; l++) 
					sum_temp += pow((tempImg[(i + k) * W + (j + l)] - l_av[(i * W + j)]), 2);
			l_value[i * W + j] = sqrt(sum_temp / pow(window_size, 2));
			g_value += l_value[i * W + j];
		}

	double thres = g_value / ((H - half_size * 2) * (W - half_size * 2));

	for (int i = half_size; i < H - half_size; i++)
		for (int j = half_size; j < W - half_size; j++) {
			int findZeroCross = 0;
			if (tempImg[(i - 1) * W + j] * tempImg[(i + 1) * W + j] < 0 && tempImg[i * W + j] == 0)
				findZeroCross = 1;
			else if (tempImg[i * W + (j - 1)] * tempImg[i * W + (j + 1)] < 0 && tempImg[i * W + j] == 0)
				findZeroCross = 1;
			else if (tempImg[i * W + j] * tempImg[(i + 1) * W + j] < 0)
				findZeroCross = 1;
			else if (tempImg[i * W + j] * tempImg[i * W + (j + 1)] > 0)
				findZeroCross = 1;

			if (findZeroCross==1)
			{
				if (l_value[i * W + j] > thres)
					Img[i * W + j] = 0;
				else
					Img[i * W + j] = 255;
			}
			else
				Img[i * W + j] = 255;
		}
}

void Q1_NonlinearLaplacian() {
	FILE* infile;

	if ((infile = fopen("Lena.pgm", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	BYTE* inImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(inImg, sizeof(BYTE), W * H, infile);
	BYTE* outImg = (BYTE*)malloc(sizeof(BYTE) * W * H);


	int* temp = (int*)malloc(sizeof(int) * W * H);
	double* l_value = (double*)malloc(sizeof(double) * W * H);

	int half_size = window_size / 2;


	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			temp[i * W + j] = 0;
			
	// 1-3
	for (int i = half_size; i < H - half_size; i++) 
		for (int j = half_size; j < W - half_size; j++) {
			BYTE max_val = 0;
			BYTE min_val = 255;
			for (int k = -half_size; k <= half_size; k++) {
				for (int l = -half_size; l <= half_size; l++) {
					max_val = max(max_val, inImg[(i + k) * W + j + l]);
					min_val = min(min_val, inImg[(i + k) * W + j + l]);
				}
			}
			temp[i * W + j] = max_val + min_val - 2 * inImg[i * W + j];
		}


	for (int i = half_size; i < H - half_size; i++)
		for (int j = half_size; j < W - half_size; j++) {
			int findZeroCross = 0;

			if (temp[(i - 1) * W + j] * temp[(i + 1) * W + j] < 0 && temp[i * W + j] == 0)
				findZeroCross = 1;
			else if (temp[i * W + (j - 1)] * temp[i * W + (j + 1)] < 0 && temp[i * W + j] == 0)
				findZeroCross = 1;
			else if (temp[i * W + j] * temp[(i + 1) * W + j] < 0)
				findZeroCross = 1;
			else if (temp[i * W + j] * temp[i * W + (j + 1)] > 0)
				findZeroCross = 1;

			if (findZeroCross == 1)
				outImg[i * W + j] = 255;
			else
				outImg[i * W + j] = 0;
		}

	// 파일 출력
	FILE* out1 = fopen("Lena_NonlinearLaplacian_1.pgm", "wb");
	fprintf(out1, "P5\n"); //magic no
	fprintf(out1, "%d %d\n", H, W);
	fprintf(out1, "%d\n", maxlevel);
	fwrite(outImg, 1, sizeof(BYTE) * H * W, out1);

	// 1-4
	Laplacian_threshold(outImg, temp, l_value);

	FILE* out2 = fopen("Lena_NonlinearLaplacian_2.pgm", "wb");
	fprintf(out2, "P5\n"); //magic no
	fprintf(out2, "%d %d\n", H, W);
	fprintf(out2, "%d\n", maxlevel);
	fwrite(outImg, 1, sizeof(BYTE) * H * W, out2);


	fclose(infile);
	fclose(out1);
	fclose(out2);
	free(inImg);
	free(outImg);
	free(temp);
	free(l_value);
}



double calc_entropy(unsigned char* image, int x, int y) {
	int hist[256] = { 0 };
	double entropy = 0.0;
	
	for (int i = -2; i <= 2; i++) {
		for (int j = -2; j <= 2; j++) {
			int idx = (y + i) * W + (x + j);
			int val = image[idx];
			hist[val]++;
		}
	}

	// 엔트로피 계산
	for (int i = 0; i < 256; i++) {
		if (hist[i] > 0) {
			double p = (double)hist[i] / 9.0;
			entropy -= p * log2(p);
		}
	}

	return entropy;
}

double calc_boundary(unsigned char* image, int x, int y) {
	// 경계값을 계산
	double boundary = 0.0;

	for (int i = -2; i <= 2; i++) {
		for (int j = -2; j <= 2; j++) {
			if (i != 0 || j != 0) {
				int idx1 = (y + i) * W + (x + j);
				int idx2 = y * W + x;
				boundary += fabs(image[idx1] - image[idx2]);
			}
		}
	}

	return boundary;
}

void Q1_EntropySketchOperator() {
	FILE* infile;
	FILE* outfile;

	if ((infile = fopen("Lena.pgm", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	BYTE* inImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	fread(inImg, sizeof(BYTE), W * H, infile);

	BYTE* outImg = (BYTE*)malloc(sizeof(BYTE) * W * H);
	memset(outImg, 0, sizeof(BYTE) * W * H);

	// 경계화소를 설정할 임계값 계산
	double max_entropy = log2(9.0);
	double max_boundary = 255.0 * 8.0;
	double threshold = (max_entropy * 0.3 + max_boundary * 0.7) / 8.0;

	// 경계화소 계산 및 출력
	for (int y = 1; y < H - 1; y++) {
		for (int x = 1; x < W - 1; x++) {
			double entropy = calc_entropy(inImg, x, y);
			double boundary = calc_boundary(inImg, x, y);
			if (entropy * 0.3 + boundary * 0.7 >= threshold) {
				outImg[y * W + x] = 0;
			}
			else {
				outImg[y * W + x] = 255;
			}
		}
	}

	outfile = fopen("Lena_EntropySketchOperator.pgm", "wb");
	fprintf(outfile, "P5\n");
	fprintf(outfile, "%d %d\n", W, H);
	fprintf(outfile, "%d\n", maxlevel);
	fwrite(outImg, sizeof(BYTE), W * H, outfile);

	fclose(outfile);
	fclose(infile);
	free(inImg);
	free(outImg);
}


void Q1_DP()
{
	FILE* infile;

	if ((infile = fopen("Lena.pgm", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	BYTE* inImg = malloc(H * W);
	BYTE* outImg = malloc(H * W);
	fread(inImg, sizeof(BYTE), H * W, infile);
	memset(outImg, 0, H * W);

	int half_size = window_size / 2;
	float* DP = malloc(H * W * sizeof(float));
	int hist[256] = { 0 };

	// DP
	int max, sum;
	for (int i = 2; i < H - 2; i++)
		for (int j = 2; j < W - 2; j++) {
			max = 0, sum = 0;
			for (int k = -2; k < 3; k++)
				for (int l = -2; l < 3; l++) {
					sum += inImg[(i + k) * W + j + l];
					if (max < inImg[(i + k) * W + j + l])
						max = inImg[(i + k) * W + j + l];
				}
			DP[i*W+j] = (float)(max - inImg[i * W + j]) / sum;
		}

	float dmax = 0, dmin = 1000;

	//Find min, max
	for (int i = 2; i < H - 2; i++)
		for (int j = 2; j < W - 2; j++) {
			if (DP[i * W + j] > dmax) dmax = DP[i * W + j];
			if (DP[i * W + j] < dmin) dmin = DP[i * W + j];
		}

	//Normalization
	for (int i = 2; i < H - 2; i++)
		for (int j = 2; j < W - 2; j++)
			outImg[i * W + j] = (BYTE)(255. * (DP[i * W + j] - dmin) / (dmax - dmin));

	//Histogram
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			hist[outImg[i * W + j]]++;

	//30%
	int cnt = 0, tmp = 255;
	for (int i = 255; i >= 0; i--) {
		cnt += hist[i];
		if (cnt > H * W * 0.3) {
			if (i != 255) tmp = i + 1;
			else tmp = 255;
			break;
		}
	}
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++) {
			if (outImg[i * W + j] >= tmp)
				outImg[i * W + j] = 255;
		}

	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++) {
			if (outImg[i * W + j] >= 255)
				outImg[i * W + j] = 0;
			else
				outImg[i * W + j] = 255;
		}

	// 파일 출력
	FILE* out = fopen("Lena_DP.pgm", "wb");
	fprintf(out, "P5\n"); //magic no
	fprintf(out, "%d %d\n", H, W);
	fprintf(out, "%d\n", maxlevel);
	fwrite(outImg, 1, sizeof(BYTE) * H * W, out);


	fclose(out);
	fclose(infile);
	free(inImg);
	free(DP);
	free(outImg);

}

void Q1_DIP()
{
	FILE* infile;

	if ((infile = fopen("Lena.pgm", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	BYTE* inImg = malloc(H * W);
	BYTE* outImg = malloc(H * W);
	fread(inImg, sizeof(BYTE), H * W, infile);
	memset(outImg, 0, H * W);

	int half_size = window_size / 2;
	float* DIP = malloc(H * W * sizeof(float));
	int hist[256] = { 0 };

	// DIP
	int max, sum;
	for (int i = 2; i < H - 2; i++)
		for (int j = 2; j < W - 2; j++) {
			max = 0, sum = 0;
			for (int k = -2; k < 3; k++)
				for (int l = -2; l < 3; l++) {
					sum += inImg[(i + k) * W + j + l];
					if (max < inImg[(i + k) * W + j + l])
						max = inImg[(i + k) * W + j + l];
				}
			DIP[i*W+j] = (float)sum / inImg[i * W + j] - (float)sum / max;
		}

	float dmax = 0, dmin = 1000;

	//Find min, max
	for (int i = 2; i < H - 2; i++)
		for (int j = 2; j < W - 2; j++) {
			if (DIP[i * W + j] > dmax) dmax = DIP[i * W + j];
			if (DIP[i * W + j] < dmin) dmin = DIP[i * W + j];
		}

	//Normalization
	for (int i = 2; i < H - 2; i++)
		for (int j = 2; j < W - 2; j++)
			outImg[i * W + j] = (BYTE)(255. * (DIP[i * W + j] - dmin) / (dmax - dmin));

	//Histogram
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			hist[outImg[i * W + j]]++;

	//30%
	int cnt = 0, tmp = 255;
	for (int i = 255; i >= 0; i--) {
		cnt += hist[i];
		if (cnt > H * W * 0.3) {
			if (i != 255) tmp = i + 1;
			else tmp = 255;
			break;
		}
	}
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++) {
			if (outImg[i * W + j] >= tmp)
				outImg[i * W + j] = 255;
		}

	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++) {
			if (outImg[i * W + j] >= 255)
				outImg[i * W + j] = 0;
			else
				outImg[i * W + j] = 255;
		}

	// 파일 출력
	FILE* out = fopen("Lena_DIP.pgm", "wb");
	fprintf(out, "P5\n"); //magic no
	fprintf(out, "%d %d\n", H, W);
	fprintf(out, "%d\n", maxlevel);
	fwrite(outImg, 1, sizeof(BYTE) * H * W, out);


	fclose(out);
	fclose(infile);
	free(inImg);
	free(outImg);

}


//Harris corner detector에 사용되는 전역 변수
float R[CH][CW];

int sx[3][3] = { -1, 0, 1, -2, 0, 2, -1, 0, 1 };
int sy[3][3] = { -1, -2, -1, 0, 0, 0, 1, 2, 1 };

int Ix2[CH][CW], Iy2[CH][CW], Ixy[CH][CW];
int Gx2[CH][CW], Gy2[CH][CW], Gxy[CH][CW];

void Q2_Harris() {
	BYTE* inImg = malloc(CW * CH);
	BYTE* outImg = malloc(CW * CH);

	FILE* infile;
	if ((infile = fopen("Ctest.raw", "rb")) == NULL) {
		printf("No Image File\n");
		return;
	}

	fread(inImg, sizeof(BYTE), CW * CH, infile);

	for (int i = 0; i < CH; i++)
		for (int j = 0; j < CW; j++)
			outImg[i * CW + j] = 0;


	// Sobel operator 적용
	int tempx, tempy;
	for (int i = 1; i < CH - 1; i++)
		for (int j = 1; j < CW - 1; j++) {
			tempx = tempy = 0;
			for (int k = -1; k < 2; k++) {
				for (int l = -1; l < 2; l++) {
					tempx += sx[k + 1][l + 1] * inImg[(i + k) * CW + j + l];
					tempy += sy[k + 1][l + 1] * inImg[(i + k) * CW + j + l];
				}
			}
			Ix2[i][j] = tempx * tempx;
			Iy2[i][j] = tempy * tempy;
			Ixy[i][j] = abs(tempx) * abs(tempy);
		}

	// 5x5 균일 가중치 평균 계산
	int sumx, sumy, sumxy;
	for (int i = 2; i < CH - 2; i++)
		for (int j = 2; j < CW - 2; j++) {
			sumx = sumy = sumxy = 0;
			for (int k = -2; k < 3; k++) {
				for (int l = -2; l < 3; l++) {
					sumx += Ix2[i + k][j + l];
					sumy += Iy2[i + k][j + l];
					sumxy += Ixy[i + k][j + l];
				}
			}
			Gx2[i][j] = sumx / 25;
			Gy2[i][j] = sumy / 25;
			Gxy[i][j] = sumxy / 25;
		}

	// Harris Corner Response Function
	for (int i = 2; i < CH - 2; i++)
		for (int j = 2; j < CW - 2; j++)
			R[i][j] = Gx2[i][j] * Gy2[i][j] - Gxy[i][j] * Gxy[i][j] - 0.05 * (Gx2[i][j] + Gy2[i][j]) * (Gx2[i][j] + Gy2[i][j]);

	// Corner 출력 (주변 인접 8pixel 비교)
	for (int i = 2; i < CH - 2; i++)
		for (int j = 2; j < CW - 2; j++) {
			if (R[i][j] > 0.01)
			{
				if (R[i][j] > R[i - 1][j - 1] && R[i][j] > R[i - 1][j] && R[i][j] > R[i - 1][j + 1] && R[i][j] > R[i][j - 1] && R[i][j] > R[i][j + 1] && R[i][j] > R[i + 1][j - 1] && R[i][j] > R[i + 1][j] && R[i][j] > R[i + 1][j + 1])
					outImg[i * CW + j] = 255;
			}
		}

	// image output
	FILE* out = fopen("Ctest_Harris.pgm", "wb");
	fprintf(out, "P5\n"); //magic no
	fprintf(out, "%d %d\n", CW, CH);
	fprintf(out, "%d\n", maxlevel);
	fwrite(outImg, sizeof(BYTE), CH * CW, out);

	fclose(infile);
	fclose(out);
	free(inImg);
	free(outImg);
}

int main()
{
	//Lena
	Q1_Sobel();
	Q1_NonlinearGradient();
	Q1_NonlinearLaplacian();
	Q1_EntropySketchOperator();
	Q1_DP();
	Q1_DIP();

	//Ctest
	Q2_Harris();
}