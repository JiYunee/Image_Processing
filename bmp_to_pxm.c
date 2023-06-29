#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <windows.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>

#define WIDTHBYTES(bits) (((bits)+31)/32*4)
#define BYTE    unsigned char
#define maxlevel 255

void BMPtoPGM()
{
	FILE* infile;

	if ((infile = fopen("ImgBW1.bmp", "rb")) == NULL) 
	{
		printf("No Image File");
		return;
	}

	// BMP Header Information
	BITMAPFILEHEADER hf;
	BITMAPINFOHEADER hInfo;
	fread(&hf, sizeof(BITMAPFILEHEADER), 1, infile); //14bytes
	fread(&hInfo, sizeof(BITMAPINFOHEADER), 1, infile);  //40 bytes

	// BMP Pallete
	RGBQUAD hRGB[256];
	fread(hRGB, sizeof(RGBQUAD), 256, infile);

	LONG W = hInfo.biWidth;
	LONG H = hInfo.biHeight;

	int rw = WIDTHBYTES(hInfo.biBitCount * W);

	// Memory
	BYTE* lpImg = malloc(hInfo.biSizeImage);
	fread(lpImg, sizeof(char), hInfo.biSizeImage, infile);
	fclose(infile);

	BYTE* pgm_RevImg = malloc(W*H);

	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			pgm_RevImg[i* W + j] = lpImg[i * rw + j];
	
	// Image Output
	FILE* outfile1 = fopen("BMPtoRevPGM.pgm", "wb");
	fprintf(outfile1, "P5\n"); //magic no
	fprintf(outfile1, "%d %d\n", W, H);
	fprintf(outfile1, "%d\n", maxlevel);
	fwrite(pgm_RevImg, sizeof(char), H * W, outfile1);

	fclose(outfile1);
	free(pgm_RevImg);

	
	BYTE* pgm_Img = malloc(W * H);

	// Image making
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			pgm_Img[(H - i - 1) * W + j] = lpImg[i * rw + j];

	// Image Output
	FILE* outfile2 = fopen("BMPtoPGM.pgm", "wb");
	fprintf(outfile2, "P5\n"); //magic no
	fprintf(outfile2, "%d %d\n", W, H);
	fprintf(outfile2, "%d\n", maxlevel);
	fwrite(pgm_Img, sizeof(char), H * W, outfile2);

	fclose(outfile2);
	free(pgm_Img);
	free(lpImg);
}

void BMPtoPPM()
{
	FILE* infile;

	if ((infile = fopen("ImgColor.bmp", "rb")) == NULL)
	{
		printf("No Image File");
		return;
	}

	// BMP Header Information
	BITMAPFILEHEADER hf;
	BITMAPINFOHEADER hInfo;
	fread(&hf, sizeof(BITMAPFILEHEADER), 1, infile); //14bytes
	fread(&hInfo, sizeof(BITMAPINFOHEADER), 1, infile);  //40 bytes

	LONG W = hInfo.biWidth;
	LONG H = hInfo.biHeight;
	int rw = WIDTHBYTES(hInfo.biBitCount * W);

	// Memory
	BYTE* lpImg = malloc(hInfo.biSizeImage);
	fread(lpImg, sizeof(char), hInfo.biSizeImage, infile);
	fclose(infile);

	BYTE* ppm_BGR = malloc(W * H * 3);

	// Image making
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			for (int k = 0; k < 3 ; k ++)
				ppm_BGR[(i*W*3) + (j*3+k)] = lpImg[(i*rw)+(j*3)+k];

	// Image Output
	FILE* outfile1 = fopen("BMPtoPPM_BGR.pgm", "wb");
	fprintf(outfile1, "P6\n"); //magic no
	fprintf(outfile1, "%d %d\n", W, H);
	fprintf(outfile1, "%d\n", maxlevel);
	fwrite(ppm_BGR, sizeof(char), H * W * 3, outfile1);

	fclose(outfile1);
	free(ppm_BGR);

	BYTE* ppm_RGB = malloc(W * H * 3);

	// Image making
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			for (int k = 0; k < 3; k++)
				ppm_RGB[(i * W * 3) + (j * 3) + (2 - k)] = lpImg[(i * rw) + (j * 3) + k];

	// Image Output
	FILE* outfile2 = fopen("BMPtoPPM_RGB.pgm", "wb");
	fprintf(outfile2, "P6\n"); //magic no
	fprintf(outfile2, "%d %d\n", W, H);
	fprintf(outfile2, "%d\n", maxlevel);
	fwrite(ppm_RGB, sizeof(char), H * W *3, outfile2);

	fclose(outfile2);
	free(ppm_RGB);

	BYTE* ppmColorImg = malloc(W * H * 3);

	// Image making
	for (int i = 0; i < H; i++)
		for (int j = 0; j < W; j++)
			for (int k = 0; k < 3; k++)
				ppmColorImg[(H-i-1)*W*3+(j*3)+(2-k)] = lpImg[(i * rw) + (j * 3) + k];

	// Image output
	FILE* outfile3 = fopen("BMPtoPPM_ColorImg.pgm", "wb");
	fprintf(outfile3, "P6\n"); //magic no
	fprintf(outfile3, "%d %d\n", W, H);
	fprintf(outfile3, "%d\n", maxlevel);
	fwrite(ppmColorImg, sizeof(char), H * W * 3, outfile3);

	fclose(outfile3);
	free(ppmColorImg);

	free(lpImg);

}

void main()
{
	BMPtoPGM();
	BMPtoPPM();
	return 0;
}
