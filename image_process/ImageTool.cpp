#include "ImageTool.h"
#include "time.h"
/**********空洞填充**********/
void ImageTool::fillHole(IplImage* srcBw, IplImage* dstBw)
{

	IplImage* Temp = cvCreateImage(cvSize(srcBw->width + 2, srcBw->height + 2), srcBw->depth, srcBw->nChannels);
	cvSet(Temp, cvScalarAll(0));
	for (int j = 1; j < Temp->height - 1; j++)
	{
		unsigned char *ptr = (unsigned char *)(srcBw->imageData + (j - 1) * srcBw->widthStep);
		unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
		for (int x = 1; x < Temp->width - 1; x++){

			ptrTemp[x] = ptr[x - 1];

		}

	}

	cvFloodFill(Temp, cvPoint(0, 0), cvScalarAll(255));
	IplImage* TempCut = cvCreateImage(cvSize(srcBw->width, srcBw->height), srcBw->depth, srcBw->nChannels);
	cvSet(TempCut, cvScalarAll(0));
	for (int j = 1; j < Temp->height - 1; j++)
	{
		unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
		unsigned char *ptrCut = (unsigned char *)(TempCut->imageData + (j - 1)* TempCut->widthStep);
		for (int x = 1; x < Temp->width - 1; x++){

			ptrCut[x - 1] = ptrTemp[x];

		}

	}
	//Mat cutImg;//裁剪延展的图像
	//Temp(Range(1, m_Size.height + 1), Range(1, m_Size.width + 1)).copyTo(cutImg);
	IplImage* TempNot = cvCreateImage(cvSize(srcBw->width, srcBw->height), srcBw->depth, srcBw->nChannels);
	cvNot(TempCut, TempNot);
	cvOr(srcBw, TempNot, dstBw);
	cvReleaseImage(&Temp);
	cvReleaseImage(&TempCut);
	cvReleaseImage(&TempNot);

}
/**********最小外接矩形，带角度**********/
CvBox2D ImageTool::smallest_rectangle2(const IplImage* In, IplImage* test, IplImage* out)
{
	//void fillHole(IplImage* srcBw, IplImage* dstBw);
	IplImage* TempTest = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
	IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
	cvSet(Temp, cvScalarAll(0));
	for (int j = 0; j < Temp->height; j++)
	{
		unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
		unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
		for (int x = 0; x < Temp->width; x++){

			ptrTemp[x] = ptr[x];

		}

	}
	CvMemStorage* storage = cvCreateMemStorage();
	CvSeq* contours = NULL;
	cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
	cvSet(Temp, cvScalarAll(0));
	cvSet(TempTest, cvScalarAll(0));
	//cvDrawContours(Temp, contours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
	CvBox2D RotatedXY = cvMinAreaRect2(contours, 0);
	
	CvPoint2D32f RotatedXYPoint[4];
	cvBoxPoints(RotatedXY, RotatedXYPoint);
	CvPoint RotatedXInt[4];
	for (int i = 0; i < 4; i++)
	{
		RotatedXInt[i].x = (int)RotatedXYPoint[i].x;
		RotatedXInt[i].y = (int)RotatedXYPoint[i].y;
	}
	for (int index = 0; index < 4; index++)
	{
		if (RotatedXInt[index].x < 0)
			RotatedXInt[index].x = 0;
		else if (RotatedXInt[index].x >In->width-2)
			RotatedXInt[index].x = In->width-1;
		if (RotatedXInt[index].y < 0)
			RotatedXInt[index].y = 0;
		else if (RotatedXInt[index].y >In->height-2)
			RotatedXInt[index].y = In->height-1;
	}
	

	cvLine(TempTest, RotatedXInt[0], RotatedXInt[1], cvScalarAll(255), 1, 8, 0);
	cvLine(TempTest, RotatedXInt[1], RotatedXInt[2], cvScalarAll(255), 1, 8, 0);
	cvLine(TempTest, RotatedXInt[2], RotatedXInt[3], cvScalarAll(255), 1, 8, 0);
	cvLine(TempTest, RotatedXInt[3], RotatedXInt[0], cvScalarAll(255), 1, 8, 0);
	fillHole(TempTest, out);
	
	cvReleaseImage(&Temp);
	cvReleaseImage(&TempTest);
	cvReleaseMemStorage(&storage);
	return RotatedXY;


}
/**********各种区域选择算子**********/
void ImageTool::select_shape(const IplImage*  In, IplImage* Out, const string &Features, float MinF, float MaxF)
{
	//////输入为黑底白色区域
	
	////根据矩形度筛选区域
	//CvBox2D smallest_rectangle2(const IplImage* In, IplImage* test, IplImage* out);
	string rectangularity = "rectangularity";
	if (!Features.compare(rectangularity)){
		////根据矩形度筛选区域

		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){

				ptrTemp[x] = ptr[x];

			}

		}
		
		//////
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int Area, AreaRectangle2;
		CvBox2D ArrayXY;
		CvMemStorage* Outstorage = cvCreateMemStorage();
		CvSeq* Outcontours = NULL;
		IplImage* Test = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		IplImage* Out1 = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Out, cvScalarAll(0));		
		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			
			cvSet(Temp, cvScalarAll(0));
			cvSet(Test, cvScalarAll(0));
			
			cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			Area = fabs(cvContourArea(IndexContours));

			ArrayXY = smallest_rectangle2(Temp, Test, Out1);

			cvFindContours(Out1, Outstorage, &Outcontours, sizeof(CvContour), CV_RETR_TREE);
			AreaRectangle2 = fabs(cvContourArea(Outcontours));

			if ((Area*1.0 / AreaRectangle2) >= MinF && (Area*1.0 / AreaRectangle2) <= MaxF)
			{
				cvDrawContours(Out, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);

			}
		}
		
		cvReleaseImage(&Temp);
		cvReleaseImage(&Test);
		cvReleaseImage(&Out1);
		cvReleaseMemStorage(&storage);
		cvReleaseMemStorage(&Outstorage);
	}
	//////////////根据面积筛选区域
	string area = "area";
	if (!Features.compare(area))
	{

		
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));		
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}

		cvSet(Out, cvScalarAll(0));
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);

		int Area;

		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			//cvSet(Temp, cvScalarAll(0));			
			//cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			Area = fabs(cvContourArea(IndexContours));

			if (Area > MinF && Area < MaxF)
			{
				cvDrawContours(Out, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);


			}


		}
		
		cvReleaseImage(&Temp);


		cvReleaseMemStorage(&storage);

	}
	///////////根据宽度筛选区域
	string width = "width";
	if (!Features.compare(width))
	{
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}
		cvSet(Out, cvScalarAll(0));
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int Area;

		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			//cvSet(Temp, cvScalarAll(0));

			//cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			//Area = fabs(cvContourArea(IndexContours));
			CvRect AreaRect= cvBoundingRect(IndexContours, 0);
			if ( AreaRect.width> MinF &&  AreaRect.width < MaxF)
			{
				cvDrawContours(Out, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);

			}


		}
		cvReleaseImage(&Temp);


		cvReleaseMemStorage(&storage);
	}
	///////////根据高度筛选区域
	string height = "height";
	if (!Features.compare(height))
	{
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}
		cvSet(Out, cvScalarAll(0));
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int Area;

		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			//cvSet(Temp, cvScalarAll(0));

			//cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			//Area = fabs(cvContourArea(IndexContours));
			CvRect AreaRect = cvBoundingRect(IndexContours, 0);
			if (AreaRect.height > MinF && AreaRect.height < MaxF)
			{
				cvDrawContours(Out, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);

			}


		}
		cvReleaseImage(&Temp);
		cvReleaseMemStorage(&storage);
	}
	string column = "column";
	if (!Features.compare(column))
	{
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int Area;
		cvSet(Out, cvScalarAll(0));
		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			//cvSet(Temp, cvScalarAll(0));

			//cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			//Area = fabs(cvContourArea(IndexContours));
			CvRect AreaRect = cvBoundingRect(IndexContours, 0);
			if ((AreaRect.width*1.0 / 2 + AreaRect.x) > MinF && (AreaRect.width*1.0 / 2 + AreaRect.x) < MaxF)
			{
				cvDrawContours(Out, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);

			}


		}
		cvReleaseImage(&Temp);


		cvReleaseMemStorage(&storage);
	}

	string row = "row";
	if (!Features.compare(row))
	{
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int Area;
		cvSet(Out, cvScalarAll(0));
		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			//cvSet(Temp, cvScalarAll(0));

			//cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			//Area = fabs(cvContourArea(IndexContours));
			CvRect AreaRect = cvBoundingRect(IndexContours, 0);
			if ((AreaRect.height*1.0 / 2 + AreaRect.y) > MinF && (AreaRect.height*1.0 / 2 + AreaRect.y) < MaxF)
			{
				cvDrawContours(Out, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);

			}


		}
		cvReleaseImage(&Temp);


		cvReleaseMemStorage(&storage);
	}


	string widthdivdeheight = "widthdivdeheight";
	if (!Features.compare(widthdivdeheight))
	{
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int Area;
		cvSet(Out, cvScalarAll(0));
		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			//cvSet(Temp, cvScalarAll(0));

			//cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			//Area = fabs(cvContourArea(IndexContours));
			CvRect AreaRect = cvBoundingRect(IndexContours, 0);
			if (AreaRect.width*1.0 / AreaRect.height > MinF && AreaRect.width*1.0 / AreaRect.height < MaxF)
			{
				cvDrawContours(Out, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);

			}


		}
		cvReleaseImage(&Temp);
		cvReleaseMemStorage(&storage);
	}






}
/**********图像相减**********/
void ImageTool::Sub(IplImage*  sub1, IplImage* sub2, IplImage* dst, int mult, int add)
{
	for (int i = 0; i < sub1->height; i++)
	{
		BYTE* ptrsub1 = (BYTE*)(sub1->imageData + i * sub1->widthStep);
		BYTE* ptrsub2 = (BYTE*)(sub2->imageData + i * sub2->widthStep);
		BYTE* ptrdst = (BYTE*)(dst->imageData + i * dst->widthStep);
		for (int j = 0; j < sub1->width; j++)
		{
			if ((ptrsub1[j] - ptrsub2[j]) * mult + add <=0){

				ptrdst[j] = 0;
			}
			else if ((ptrsub1[j] - ptrsub2[j]) * mult + add >=255){
				ptrdst[j] = 255;
			}
			else{
				ptrdst[j] = (ptrsub1[j] - ptrsub2[j]) * mult + add;
			}


		}
	}

}
/**********傅立叶变换**********/
void ImageTool::IDFT(IplImage * fourier, IplImage* dst) {
	//IplImage* dst = cvCreateImage(cvGetSize(fourier), IPL_DEPTH_8U, 1);
	int dft_H, dft_W;
	dft_H = fourier->height;
	dft_W = fourier->width;
	CvMat *dst_Re = cvCreateMat(dft_H, dft_W, CV_64FC1);
	// double Re, Im;     
	CvMat *dst_Im = cvCreateMat(dft_H, dft_W, CV_64FC1);
	//Imaginary part      
	CvMat *sum_dst = cvCreateMat(dft_H, dft_W, CV_64FC2);
	//2 channels (dst_Re, dst_Im)  
	CvMat *sum_src = cvCreateMat(dft_H, dft_W, CV_64FC2);
	cvConvert(fourier, sum_src);
	cvDFT(sum_src, sum_dst, CV_DXT_INV_SCALE, 0);
	cvSplit(sum_dst, dst_Re, dst_Im, 0, 0);
	cvConvert(dst_Re, dst);
	cvReleaseMat(&dst_Re);
	cvReleaseMat(&dst_Im);
	cvReleaseMat(&sum_src);
	cvReleaseMat(&sum_dst);
	//return dst;
}
void ImageTool::BuildDFTImage(IplImage *fourier, IplImage *dst) {
	IplImage *image_Re = 0, *image_Im = 0;
	image_Re = cvCreateImage(cvGetSize(fourier), IPL_DEPTH_64F, 1);
	image_Im = cvCreateImage(cvGetSize(fourier), IPL_DEPTH_64F, 1);
	//Imaginary part   
	cvSplit(fourier, image_Re, image_Im, 0, 0);
	// Compute the magnitude of the spectrum Mag = sqrt(Re^2 + Im^2)   
	cvPow(image_Re, image_Re, 2.0);
	cvPow(image_Im, image_Im, 2.0);
	cvAdd(image_Re, image_Im, image_Re);
	cvPow(image_Re, image_Re, 0.5);
	cvReleaseImage(&image_Im);
	cvAddS(image_Re, cvScalar(1.0), image_Re);
	// 1 + Mag  
	cvLog(image_Re, image_Re);// log(1 + Mag)    
	//重新安排傅里叶图像中心 
	// Rearrange the quadrants of Fourier image so that the origin is at
	// the image center 
	double minVal = 0, maxVal = 0;
	cvMinMaxLoc(image_Re, &minVal, &maxVal); // Localize minimum and maximum values   
	CvScalar min;
	min.val[0] = minVal;
	double scale = 255 / (maxVal - minVal);
	cvSubS(image_Re, min, image_Re);
	cvConvertScale(image_Re, dst, scale);
	cvReleaseImage(&image_Re);
	// Rearrange the quadrants of Fourier image so that the origin is at
	// the image center  
	int nRow, nCol, i, j, cy, cx;
	uchar tmp13, tmp24;
	nRow = fourier->height;
	nCol = fourier->width;
	cy = nRow / 2; // image center  
	cx = nCol / 2;
	for (j = 0; j < cy; j++)  {
		for (i = 0; i < cx; i++)
		{
			tmp13 = CV_IMAGE_ELEM(dst, uchar, j, i);
			CV_IMAGE_ELEM(dst, uchar, j, i) = CV_IMAGE_ELEM(dst, uchar, j + cy, i + cx);
			CV_IMAGE_ELEM(dst, uchar, j + cy, i + cx) = tmp13;
			tmp24 = CV_IMAGE_ELEM(dst, uchar, j, i + cx);
			CV_IMAGE_ELEM(dst, uchar, j, i + cx) = CV_IMAGE_ELEM(dst, uchar, j + cy, i);
			CV_IMAGE_ELEM(dst, uchar, j + cy, i) = tmp24;
		}
	}
}

void ImageTool::PassFilter(IplImage * fourier, int FLAG, double D0, int n1) {
	int i, j;
	int state = -1;
	double tempD;
	long width, height;
	width = fourier->width;
	height = fourier->height;
	long x, y;  x = width / 2;
	y = height / 2;
	CvMat* H_mat;
	H_mat = cvCreateMat(fourier->height, fourier->width, CV_64FC2);
	for (i = 0; i < height; i++){
		for (j = 0; j < width; j++){
			if (i > y && j > x){
				state = 3;
			}
			else if (i > y){
				state = 1;
			}
			else if (j > x){
				state = 2;
			}
			else{
				state = 0;
			}
			switch (state){
			case 0:      tempD = (double)sqrt(1.0*i * i + j * j);
				break;
			case 1:      tempD = (double)sqrt(1.0*(height - i) * (height - i) + j * j);
				break;
			case 2:      tempD = (double)sqrt(1.0*i * i + (width - j) * (width - j));
				break;
			case 3:      tempD = (double)sqrt(1.0*(height - i) * (height - i) + (width - j) * (width - j));
				break;
			default:
				break;
			}
			switch (FLAG){
			case 0:
				if (tempD <= D0){
					((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2] = 1.0;
					((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2 + 1] = 0.0;
				}
				else{
					((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2] = 0.0;
					((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2 + 1] = 0.0;
				}
				break;
			case 1:
				if (tempD <= D0){
					((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2] = 0.0;
					((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2 + 1] = 0.0;
				}
				else{
					((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2] = 1.0;
					((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2 + 1] = 0.0;
				}
				break;
			case 2:
				tempD = 1 / (1 + pow(tempD / D0, 2 * n1));
				((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2] = tempD;
				((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2 + 1] = 0.0;
				break;
			case 3:
				tempD = 1 / (1 + pow(D0 / tempD, 2 * n1));
				((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2] = tempD;
				((double*)(H_mat->data.ptr + H_mat->step * i))[j * 2 + 1] = 0.0;
				break;
			default:
				break;
			}
		}
	}
	cvMulSpectrums(fourier, H_mat, fourier, CV_DXT_ROWS);
	cvReleaseMat(&H_mat);
}
void  ImageTool::DFT(IplImage * src, IplImage * fourier) {
	//IplImage* fourier = cvCreateImage(cvGetSize(src), IPL_DEPTH_64F, 2);
	int dft_H, dft_W;
	dft_H = src->height;
	dft_W = src->width;
	CvMat *src_Re = cvCreateMat(dft_H, dft_W, CV_64FC1);
	// double Re, Im;     
	CvMat *src_Im = cvCreateMat(dft_H, dft_W, CV_64FC1);
	//Imaginary part  
	CvMat *sum_src = cvCreateMat(dft_H, dft_W, CV_64FC2);
	//2 channels (src_Re, src_Im)      
	CvMat *sum_dst = cvCreateMat(dft_H, dft_W, CV_64FC2);
	//2 channels (dst_Re, dst_Im)  
	cvConvert(src, src_Re);
	cvZero(src_Im);
	cvMerge(src_Re, src_Im, 0, 0, sum_src);
	cvDFT(sum_src, sum_dst, CV_DXT_FORWARD, 0);
	cvConvert(sum_dst, fourier);
	cvReleaseMat(&src_Re);
	cvReleaseMat(&src_Im);
	cvReleaseMat(&sum_src);
	cvReleaseMat(&sum_dst);
	//return fourier;
}
/*******************************************************************/
/**********图像复制**********/
void  ImageTool::CloneImage(IplImage * In, IplImage * Out){

	if (In->nChannels == 1)
	{

		cvSet(Out, cvScalarAll(0));
		for (int j = 0; j < Out->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Out->imageData + (j)* Out->widthStep);
			for (int x = 0; x < Out->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}
	}
	else if (In->nChannels == 3)
	{
		cvSet(Out, cvScalarAll(0));
		for (int j = 0; j < Out->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Out->imageData + (j)* Out->widthStep);
			for (int x = 0; x < Out->width; x++){
				ptrTemp[3 * x + 1] = ptr[3 * x + 1];
				ptrTemp[3 * x + 2] = ptr[3 * x + 2];
				ptrTemp[3 * x + 3] = ptr[3 * x + 3];
			}

		}
	}

}
int  ImageTool::SumArray(int *a,int length){
	int sum = 0;
	for (int i = 0; i < length; i++){
		sum = sum + a[i];

	}
	return sum;


}
/**********获取文件夹文件名**********/
vector<string> ImageTool::get_all_files_names_within_folder(string folder)
{
	vector<string> names;
	WIN32_FIND_DATA fd;
	HANDLE hFind = ::FindFirstFile(folder.append("\\*").c_str(), &fd);
	if (hFind != INVALID_HANDLE_VALUE) {
		do {

			if (fd.cFileName[0] == '.' || fd.cFileName[0] == '..')
				continue;

			if ((fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) || fd.cFileName[0] == '.') {
				names.push_back(fd.cFileName);
			}
		} while (::FindNextFile(hFind, &fd));
		::FindClose(hFind);
	}
	return names;
}
/**********获取图片**********/
void ImageTool::getFiles(string path, vector<string>& files, string imagetype)
{
	//文件句柄  
	long   hFile = 0;
	//文件信息  
	struct _finddata_t fileinfo;
	string p;
	if ((hFile = _findfirst(p.assign(path).append("\\*." + imagetype).c_str(), &fileinfo)) != -1)
	{
		do
		{
			if ((fileinfo.attrib &  _A_SUBDIR))
			{
				if (strcmp(fileinfo.name, ".") != 0 && strcmp(fileinfo.name, "..") != 0)
					getFiles(p.assign(path).append("\\").append(fileinfo.name), files, imagetype);
			}
			else
			{
				files.push_back(p.assign(path).append("\\").append(fileinfo.name));
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}

/**********截取小图**********/
void ImageTool::CropImage(IplImage * in, IplImage * out, int row_left, int col_left, int row_right, int col_right,int channels)
{
	if (channels == 3)
	{

		int testy = 0;
		int testx = 0;
		for (int i = row_left; i < row_right; i++)
		{
			//BYTE* ptr = (BYTE*)(InImageTemp1->imageData + i*InImageTemp1->widthStep);
			BYTE* ptrin = (BYTE*)(in->imageData + i*in->widthStep);
			BYTE* ptrout = (BYTE*)(out->imageData + testx*out->widthStep);
			testy = 0;
			for (int j = col_left; j < col_right; j++)
			{
				ptrout[3 * testy] = ptrin[3 * j];
				ptrout[3 * testy + 1] = ptrin[3 * j + 1];
				ptrout[3 * testy + 2] = ptrin[3 * j + 2];


				testy++;
			}
			testx++;
		}
	}
	if (channels == 1)
	{

		int testy = 0;
		int testx = 0;
		for (int i = row_left; i < row_right; i++)
		{
			//BYTE* ptr = (BYTE*)(InImageTemp1->imageData + i*InImageTemp1->widthStep);
			BYTE* ptrin = (BYTE*)(in->imageData + i*in->widthStep);
			BYTE* ptrout = (BYTE*)(out->imageData + testx*out->widthStep);
			testy = 0;
			for (int j = col_left; j < col_right; j++)
			{
				ptrout[testy] = ptrin[j];
				testy++;
			}
			testx++;
		}
	}

}

/**********图像增强**********/
void ImageTool::scale_image(IplImage * in, IplImage * out, int Mult, int Add)
{
	float LowerLimit = 0.0;
	float UpperLimit = 255.0;
	if (in->nChannels == 3)
	{

		for (int i = 0; i < in->height; i++)
		{

			BYTE* ptrin = (BYTE*)(in->imageData + i*in->widthStep);
			BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);

			for (int j = 0; j < in->width; j++)
			{
				if (ptrin[3 * j] * Mult + Add>255)
					ptrout[3 * j] = 255;
				else
					ptrout[3 * j] = ptrin[3 * j] * Mult + Add;
				    				
				if (ptrin[3 * j+1] * Mult+Add > 255)
					ptrout[3 * j+1] = 255;
				else
					ptrout[3 * j+1] = ptrin[3 * j+1] * Mult + Add;

				if (ptrin[3 * j + 2] * Mult + Add > 255)
					ptrout[3 * j + 2] = 255;
				else
					ptrout[3 * j + 2] = ptrin[3 * j + 2] * Mult + Add;

			}
		}
	}
	else if (in->nChannels == 1)
	{

		for (int i = 0; i < in->height; i++)
		{

			BYTE* ptrin = (BYTE*)(in->imageData + i*in->widthStep);
			BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
			for (int j = 0; j < in->width; j++)
			{
				ptrout[j] = ptrin[j] * Mult + Add;

			}
		}


	}
	//Calculate scaling parameters
	//float Mult = (UpperLimit - LowerLimit) / (MaxGray - MinGray);
	//float Add = -Mult * MinGray + LowerLimit;
	
}
/**********动态二值化**********/
void ImageTool::Dyn_Threhold(IplImage * in, IplImage * mean, IplImage * out, int offset, int flag)
{
	//// 0 dark 1 light
	if (flag == 0)
	{

		for (int i = 0; i < in->height; i++)
		{

			BYTE* ptrin = (BYTE*)(in->imageData + i*in->widthStep);
			BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
			BYTE* ptrmean = (BYTE*)(mean->imageData + i*mean->widthStep);
			for (int j = 0; j < in->width; j++)
			{
				if (ptrin[j] <= ptrmean[j] - offset)
					ptrout[j] = 255;
				else
					ptrout[j] = 0;
			}
		}
	}
	if (flag == 1)
	{

		for (int i = 0; i < in->height; i++)
		{

			BYTE* ptrin = (BYTE*)(in->imageData + i*in->widthStep);
			BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
			BYTE* ptrmean = (BYTE*)(mean->imageData + i*mean->widthStep);
			for (int j = 0; j < in->width; j++)
			{
				if (ptrin[j] >= ptrmean[j] + offset)
					ptrout[j] = 255;
				else
					ptrout[j] = 0;
			}
		}
	}



}
/**********截取外接矩形，不带角度**********/
void ImageTool::Gen_Rectangle1(IplImage * out, int Row1, int Column1, int Row2, int Column2)
{
	if (out->nChannels == 3)
	{
		//cvSet(out, cvScalarAll(0));
		for (int i = Row1; i <=Row2; i++)
		{
			//BYTE* ptrin = (BYTE*)(out->imageData + i*out->widthStep);
			if (i == Row2)
			{
				BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
				for (int j = Column1; j <=Column2; j++)
				{
						ptrout[3 * j] = 255;
						ptrout[3 * j + 1] = 255;
						ptrout[3 * j + 2] = 255;

				}
			}
			if (i == Row1)
			{
				BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
				for (int j = Column1; j <= Column2; j++)
				{

						ptrout[3 * j] = 255;
						ptrout[3 * j + 1] = 255;
						ptrout[3 * j + 2] = 255;

				}
			}
		
			if (i > Row1 && i<Row2)
			{
				BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
				for (int j = Column1; j <= Column2; j++)
				{
					if (j == Column1 || j == Column2)
					{
						ptrout[3 * j] = 255;
						ptrout[3 * j + 1] = 255;
						ptrout[3 * j + 2] = 255;
					}

				}
			}	
		}		
	}
	else if(out->nChannels == 1)
	{
		//cvSet(out, cvScalarAll(0));
		for (int i = Row1; i <= Row2; i++)
		{
			//BYTE* ptrin = (BYTE*)(out->imageData + i*out->widthStep);
			if (i == Row2)
			{
				BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
				for (int j = Column1; j <= Column2; j++)
				{
					ptrout[j] = 255;
				}
			}
			if (i == Row1)
			{
				BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
				for (int j = Column1; j <= Column2; j++)
				{
					ptrout[j] = 255;
				}
			}

			if (i > Row1 && i < Row2)
			{
				BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
				for (int j = Column1; j <= Column2; j++)
				{
					if (j == Column1 || j == Column2)
					{
						ptrout[j] = 255;
					}

				}
			}



		}
	}



}
/**********获取最大区域**********/
int ImageTool::FindMaxArea(IplImage * in, IplImage * out)
{
	CvMemStorage *storage = cvCreateMemStorage();
	CvSeq *contours = NULL;
	int MaxArea;
	MaxArea = 0;
	IplImage *inTemp = cvCreateImage(cvSize(in->width, in->height), in->depth, in->nChannels);
	CloneImage(in,inTemp);
	cvFindContours(inTemp, storage, &contours, sizeof(CvContour), CV_RETR_TREE, CV_CHAIN_APPROX_NONE);
	cvSet(out, cvScalarAll(0));
	for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
	{
		if (fabs(cvContourArea(IndexContours)) > MaxArea)
		{

			MaxArea = fabs(cvContourArea(IndexContours));
			cvSet(out, cvScalarAll(0));
			cvDrawContours(out, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
		}
	}
	cvReleaseMemStorage(&storage);
	cvReleaseImage(&inTemp);
	return MaxArea;
}
/**********获取区域特征**********/
float ImageTool::RegionFeatures(IplImage * In, const string &Features)
{
	string area = "area";
	if (!Features.compare(area)){
		CvMemStorage *storage = cvCreateMemStorage();
		CvSeq *contours = NULL;
		int Area;
		Area = 0;
		IplImage *inTemp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		CloneImage(In, inTemp);
		cvFindContours(inTemp, storage, &contours, sizeof(CvContour), CV_RETR_TREE, CV_CHAIN_APPROX_NONE);
		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			Area = Area+fabs(cvContourArea(IndexContours));
		}
		cvReleaseMemStorage(&storage);
		cvReleaseImage(&inTemp);
		return Area;
	}
	string width = "width";
	if (!Features.compare(width))
	{
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		//cvSet(Out, cvScalarAll(0));
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int Width=0;

		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			//cvSet(Temp, cvScalarAll(0));

			//cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			//Area = fabs(cvContourArea(IndexContours));
			CvRect AreaRect = cvBoundingRect(IndexContours, 0);
			Width = Width + AreaRect.width;


		}
		
		cvReleaseImage(&Temp);


		cvReleaseMemStorage(&storage);
		return Width;
	}
	///////////根据高度筛选区域
	string height = "height";
	if (!Features.compare(height))
	{
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		//cvSet(Out, cvScalarAll(0));
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){
				ptrTemp[x] = ptr[x];

			}

		}
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int Height=0;
		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			//cvSet(Temp, cvScalarAll(0));

			//cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			//Area = fabs(cvContourArea(IndexContours));
			CvRect AreaRect = cvBoundingRect(IndexContours, 0);
			Height = Height + AreaRect.height;
		}
		
		cvReleaseImage(&Temp);
		cvReleaseMemStorage(&storage);
		return Height;
	}
	string rectangularity = "rectangularity";
	
	if (!Features.compare(rectangularity)){
		////根据矩形度筛选区域
		CvMemStorage* storage = cvCreateMemStorage();
		CvSeq* contours = NULL;
		/////复制图片，给findcontour使用
		IplImage* Temp = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		cvSet(Temp, cvScalarAll(0));
		for (int j = 0; j < Temp->height; j++)
		{
			unsigned char *ptr = (unsigned char *)(In->imageData + (j)* In->widthStep);
			unsigned char *ptrTemp = (unsigned char *)(Temp->imageData + (j)* Temp->widthStep);
			for (int x = 0; x < Temp->width; x++){

				ptrTemp[x] = ptr[x];

			}

		}
		//////
		cvFindContours(Temp, storage, &contours, sizeof(CvContour), CV_RETR_TREE);
		int coutrec = 0;
		float rectangularity = 0;
		int Area = 0, AreaRectangle2 = 0;
		CvBox2D ArrayXY;
		CvMemStorage* Outstorage = cvCreateMemStorage();
		CvSeq* Outcontours = NULL;
		IplImage* Test = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		IplImage* Out1 = cvCreateImage(cvSize(In->width, In->height), In->depth, In->nChannels);
		
		for (CvSeq* IndexContours = contours; IndexContours != NULL; IndexContours = IndexContours->h_next)
		{
			cvSet(Temp, cvScalarAll(0));
			cvSet(Test, cvScalarAll(0));
			cvDrawContours(Temp, IndexContours, cvScalar(255, 255, 255), cvScalar(255, 255, 255), 0, CV_FILLED, 8);
			Area = fabs(cvContourArea(IndexContours));
			ArrayXY = smallest_rectangle2(Temp, Test, Out1);

			cvFindContours(Out1, Outstorage, &Outcontours, sizeof(CvContour), CV_RETR_TREE, CV_CHAIN_APPROX_NONE);
			AreaRectangle2 = fabs(cvContourArea(Outcontours));
			
			rectangularity = rectangularity+Area*1.0 / AreaRectangle2;
			coutrec++;
		}
		rectangularity = rectangularity*1.0 / coutrec;
		cvReleaseImage(&Temp);
		cvReleaseImage(&Test);
		cvReleaseImage(&Out1);
		cvReleaseMemStorage(&storage);
		cvReleaseMemStorage(&Outstorage);
		return rectangularity;
	}
	
}
/**********获取两区域不同区域**********/
void ImageTool::RegionDifference(IplImage * in1, IplImage * in2, IplImage * out)
{
	cvSet(out,cvScalarAll(0));
	for (int i = 0; i < in1->height; i++)
	{

		BYTE* ptrin1 = (BYTE*)(in1->imageData + i*in1->widthStep);
		BYTE* ptrout = (BYTE*)(out->imageData + i*out->widthStep);
		BYTE* ptrin2 = (BYTE*)(in2->imageData + i*in2->widthStep);
		for (int j = 0; j < in1->width; j++)
		{
			//if (ptrin1[j] != ptrin2[j] && (ptrin1[j] == 255 || ptrin2[j]==255))
			if (ptrin1[j] != ptrin2[j] && ptrin1[j] == 255)
				ptrout[j] = 255;
			else
				ptrout[j] = 0;
		}
	}



}
ImageTool::~ImageTool()
{
	
}
ImageTool::ImageTool()
{
	
}

