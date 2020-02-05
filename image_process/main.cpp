#include "ImageTool.h"
#include "string"
#include <vector> ; 
#include "time.h"
using namespace std;
int main(int argc, char** argv){
	
	string FileNamegeng = "D:\\Image\\";
	ImageTool ImageToolClass;
	ImgQuImprove ImgQuImproveClass;
	vector<string> Filesgeng;

	Filesgeng = ImageToolClass.get_all_files_names_within_folder(FileNamegeng);
	for (int indexfiles = 0; indexfiles < Filesgeng.size(); indexfiles++)
	{
		string ImagePath = FileNamegeng + "\\" + Filesgeng[indexfiles];
		vector<string>ImageNames;
		string imagetype = "jpg";
		ImageToolClass.getFiles(ImagePath, ImageNames, imagetype);
		for (int IndexImages =0; IndexImages < ImageNames.size(); IndexImages++)
		{
			IplImage *OrgImage = cvLoadImage(ImageNames[IndexImages].c_str(), 1);
			IplImage* ImproveImg = cvCreateImage(cvSize(OrgImage->width, OrgImage->height), OrgImage->depth, OrgImage->nChannels);



			cvReleaseImage(&ImproveImg);
			cvReleaseImage(&OrgImage);
		}
	}












}