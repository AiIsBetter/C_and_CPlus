#include "highgui.h"
#include "cv.h"
#include "cxcore.h"
#include <iosfwd>  
#include <memory>  
#include <utility>  
#include <vector>  
#include <iostream>  
#include <string.h>  
#include <sstream>  
#include <io.h>
#include <windows.h>
#include <tchar.h>
#include "cxcore.h"
using namespace std;

class ImageTool
{
public:
	ImageTool();
	virtual ~ImageTool();
	void fillHole(IplImage* srcBw, IplImage* dstBw);
	CvBox2D smallest_rectangle2(const IplImage* In, IplImage* test, IplImage* out);
	void select_shape(const IplImage*  In, IplImage* Out, const string &Features, float MinF, float MaxF);
	

	void Sub(IplImage*  sub1, IplImage* sub2, IplImage* dst, int mult, int add);

	void BuildDFTImage(IplImage *fourier, IplImage *dst);
	void PassFilter(IplImage * fourier, int FLAG, double D0, int n1);
	void DFT(IplImage * src, IplImage * fourier);
	void IDFT(IplImage * fourier, IplImage* dst);
	void CloneImage(IplImage * In, IplImage * Out);
	int  SumArray(int *a, int length);
	vector<string> get_all_files_names_within_folder(string folder);
	void getFiles(string path, vector<string>& files, string imagetype);
	void CropImage(IplImage * in, IplImage * out, int row_left, int col_left, int row_right, int col_right, int channels);
	void scale_image(IplImage * in, IplImage * out, int Mult, int Add);
	void Dyn_Threhold(IplImage * in, IplImage * mean, IplImage * out, int offset, int flag);
	void Gen_Rectangle1(IplImage * out, int Row1, int Column1, int Row2, int Column2);
	int FindMaxArea(IplImage * in, IplImage * out);
	float RegionFeatures(IplImage * in, const string &Features);
	void RegionDifference(IplImage * in1, IplImage * in2, IplImage * out);
};