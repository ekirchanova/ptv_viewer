#include "image.h"

#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include <opencv2/opencv.hpp>
#include <math.h>


int prepareCalibImageCV(const char * fileName, camera_sequence *c){
    using namespace cv;
    using namespace std;
    String str(fileName);
    Mat  image= imread(str,IMREAD_GRAYSCALE);


  printf("reading image %s w=%d h=%d\n",fileName,image.size().width,image.size().height);
    Size boardSize(19,19);
    vector<Point2f> pointbuf;
    findCirclesGrid( image, boardSize, pointbuf );
    c->imagePoints.push_back(pointbuf);
    printf("found %d points\n",pointbuf.size());

    double dr = 0.051;
    int ni = 19;
    vector<Point3f> singleObj;
    for( int i = 0; i < ni; i++ )
        for( int j = 0; j < ni; j++ )
            singleObj.push_back(Point3f(float(dr*(j - (ni-1)*0.5)),
                                        float(dr*(i - (ni-1)*0.5)), 0));
    c->objectPoints.push_back(singleObj);


    c->images.push_back(image);
    /*for (int i=0; i<pointbuf.size(); i++){
        printf("point %d x=%f y=%f \n",i,pointbuf[i].x,pointbuf[i].y);
        c->centers[i][0]=pointbuf[i].x*1.0/c->w;
        c->centers[i][1]=1.0 - pointbuf[i].y*1.0/c->h;
        c->centers[i][2]=0.0;
    }
    c->centers_size = pointbuf.size();

    /*
    //filling c TODO pixel format conversion
    c->h = image.size().height;
    c->w = image.size().width;
    int totalSize = c->h*c->w*3;
    c->data=(byte*)realloc((void*)(c->data), totalSize);
    printf("realloc done h=%d w=%d \n",c->h,c->w);
    //  memcpy((void*)(c.data),(void*)(pixels),totalSize);


    int index=0;
    for(int j = 0; j < c->h ; ++j)
        for(int i = 0; i < c->w ; ++i){

            unsigned char B,R,G;
            B = image.at<unsigned char>(c->h-1-j,i);
            G = image.at<unsigned char>(c->h-1-j,i);
            R = image.at<unsigned char>(c->h-1-j,i);

            c->data[index] = R;
            c->data[index+1] = G;
            c->data[index+2] = B;
            index+=3;
        }*/
    //contours



    /*
    Mat res(image.size(), CV_8UC1);
    Canny(image, res, 100 , 200, 3);

    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
     findContours(res, contours, hierarchy, RETR_TREE, CHAIN_APPROX_SIMPLE); // Find contours with hierarchy

    vector<Moments> mu_tbl(contours.size());
    vector<Point2f> mc_tbl(contours.size());
    int num=0;
    for( int i = 0; i < contours.size(); i++ )
    {	mu_tbl[i] = moments( contours[i], false );
        //mc_tbl[i] = Point2f( mu_tbl[i].m10/mu_tbl[i].m00 , mu_tbl[i].m01/mu_tbl[i].m00 );
        //  printf("cont %d: rad=%f size=%d\n",i, mu_tbl[i].m00,contours[i].size());
        for( int j = 0; j < contours[i].size(); j++ ){
            //  printf("%d %f %f \n",num,contours[i][j].x*1.0/c->w,1.0 - contours[i][j].y*1.0/c->h);
            c->points[num][0] = contours[i][j].x*1.0/c->w;
            c->points[num][1] = 1.0 - contours[i][j].y*1.0/c->h;
            c->points[num][2] = 0.0;
            if(num <100000)num++;
        }
    }
    c->pointSize = num;*/
    return 1;
}

void try_calibration(camera_sequence *c){
    using namespace cv;
    using namespace std;
    c->distCoeffs = Mat::zeros(8, 1, CV_64F);
    double rms = calibrateCamera(c->objectPoints, c->imagePoints, c->images[0].size(), c->cameraMatrix,
                                 c->distCoeffs, c->rvects, c->tvects,  CALIB_USE_LU);
    printf("%d rvecs has %d rows and %d cols \n",c->rvects.size(),c->rvects[0].rows, c->rvects[0].cols);
    printf("%d tvecs has %d rows and %d cols \n",c->tvects.size(),c->tvects[0].rows, c->tvects[0].cols);
    printf("RMS error reported by calibrateCamera: %g\n %d %d", rms,c->cameraMatrix.rows,c->cameraMatrix.cols);
    printf("mat:\n %f %f %f \n %f %f %f \n %f %f %f\n",c->cameraMatrix.at<double>(0,0),c->cameraMatrix.at<double>(0,1),c->cameraMatrix.at<double>(0,2),
           c->cameraMatrix.at<double>(1,0),c->cameraMatrix.at<double>(1,1),c->cameraMatrix.at<double>(1,2),
           c->cameraMatrix.at<double>(2,0),c->cameraMatrix.at<double>(2,1),c->cameraMatrix.at<double>(2,2));
   /*
    c_[0].rvects.clear();
    c_[0].tvects.clear();

    for (int i=0; i<rvecs.size(); ++i){
        c_[0].rvects.push_back(rvecs[i]);
        c_[0].tvects.push_back(tvecs[i]);
    }
    c_[0].cameraMatrix = cameraMatrix;
    c_[0].distCoeffs = distCoeffs;*/


    OpenGLMatrix modelv;
    convertOpenCVProjectionMatrixToOpenGL(c->cameraMatrix,0.1,1.5,c->projM.m);
    for (int i=0;i<c->rvects.size(); ++i){
    convertOpenCVModelViewMatrixToOpenGL(c->rvects[i],c->tvects[i],modelv.m);
    c->modelvMatrices.push_back(modelv);
    }

}

void convertOpenCVModelViewMatrixToOpenGL(cv::Mat rvec, cv::Mat tvec, GLdouble* GLModelV ){
    cv::Mat rotation;
    cv::Mat viewMatrix = cv::Mat::zeros(4, 4, CV_64FC1);
    cv::Rodrigues(rvec, rotation);

    for(unsigned int row=0; row<3; ++row)
    {
        for(unsigned int col=0; col<3; ++col)
        {
            viewMatrix.at<double>(row, col) = rotation.at<double>(row, col);
        }
        viewMatrix.at<double>(row, 3) = tvec.at<double>(row, 0);
    }
    viewMatrix.at<double>(3, 3) = 1.0f;

    cv::Mat cvToGl = cv::Mat::zeros(4, 4, CV_64F);
    cvToGl.at<double>(0, 0) = 1.0f;
    cvToGl.at<double>(1, 1) = -1.0f; // Invert the y axis
    cvToGl.at<double>(2, 2) = -1.0f; // invert the z axis
    cvToGl.at<double>(3, 3) = 1.0f;
    viewMatrix = cvToGl * viewMatrix;

    cv::Mat glViewMatrix = cv::Mat::zeros(4, 4, CV_64F);
    cv::transpose(viewMatrix , glViewMatrix);
    //glMatrixMode(GL_MODELVIEW);
    //glLoadMatrixd(&glViewMatrix.at<double>(0, 0));
    int n=0;
    for (int i=0;i<4; ++i)
        for (int j=0; j<4; ++j)
        {
            GLModelV[n] = glViewMatrix.at<double>(i,j);
            n++;
        }
}


void convertOpenCVProjectionMatrixToOpenGL(cv::Mat cameraMatrix, double zN,double zF,GLdouble* GLProjM){

    double fx,fy,W_,H_;
    fx = cameraMatrix.at<double>(0,0);
    fy = cameraMatrix.at<double>(1,1);
    W_ = cameraMatrix.at<double>(0,2)*2.0;
    H_ = cameraMatrix.at<double>(1,2)*2.0;
    GLProjM[0] = 2*fx/W_; GLProjM[4] = 0.0;     GLProjM[8] = 0.0; GLProjM[12] = 0.0;
    GLProjM[1] = 0.0;     GLProjM[5] = 2*fy/H_; GLProjM[9] = 0.0; GLProjM[13] = 0.0;
    GLProjM[2] = 0.0;     GLProjM[6] = 0.0;     GLProjM[10] = -(zF+zN)/(zF-zN); GLProjM[14] = -2.0*(zF*zN)/(zF-zN);
    GLProjM[3] = 0.0;     GLProjM[7] = 0.0;     GLProjM[11] = -1.0; GLProjM[15] = 0.0;
}

void undistortPoints(size_t numberImage, camera_sequence *c)
{
    using namespace cv;
    using namespace std;

    vector<Point2f> pointbuf;
    undistortPoints(c->imagePoints[numberImage], pointbuf, c->cameraMatrix, c->distCoeffs);
    c->undistortImagePoints.push_back(pointbuf);
    printf("found %d points\n",pointbuf.size());
}

/*
int loadBmpImage( const char * fileName, camera_sequence *c)
{
        //Open the file for reading in binary mode
        FILE *imageFile = fopen(fileName, "rb");
        if ( imageFile == NULL ){
         printf("no such file! %s \n",fileName);
            return 0;
        }
        //Read data offset
        int32 dataOffset;
        fseek(imageFile, DATA_OFFSET_OFFSET, SEEK_SET);
        fread(&dataOffset, 4, 1, imageFile);
        //Read width
        fseek(imageFile, WIDTH_OFFSET, SEEK_SET);
        int32 width,height;
        fread(&width, 4, 1, imageFile);
        //Read height
        fseek(imageFile, HEIGHT_OFFSET, SEEK_SET);
        fread(&height, 4, 1, imageFile);
        //Read bits per pixel
        int16 bitsPerPixel;
        fseek(imageFile, BITS_PER_PIXEL_OFFSET, SEEK_SET);
        fread(&bitsPerPixel, 2, 1, imageFile);
        //Allocate a pixel array
        int32 bytesPerPixel;
        bytesPerPixel = ((int32)bitsPerPixel) / 8;


        //Rows are stored bottom-up
        //Each row is padded to be a multiple of 4 bytes.
        //We calculate the padded row size in bytes
        int paddedRowSize = (int)(4 * ceil((float)(width) / 4.0f))*(bytesPerPixel);
        //We are not interested in the padded bytes, so we allocate memory just for
        //the pixel data
        int unpaddedRowSize = (width)*(bytesPerPixel);
        //Total size of the pixel data in bytes
        int totalSize = unpaddedRowSize*(height);
        byte* pixels;
        pixels = (byte*)malloc(totalSize);
        //Read the pixel data Row by Row.
        //Data is padded and stored bottom-up
        int i = 0;
        //point to the last row of our pixel array (unpadded)
        printf("image width = %d image height = %d bpp=%d bypp=%d prs=%d nprs=%d totalSize=%d \n",width,height,bitsPerPixel, bytesPerPixel, paddedRowSize, unpaddedRowSize,totalSize);
        fseek(imageFile, dataOffset, SEEK_SET);
        fread((void*)pixels,totalSize,1,imageFile);


        fclose(imageFile);
        printf("file read \n");
        //filling c TODO pixel format conversion
        c->h = height;
        c->w = width;
        c->data=(byte*)realloc((void*)(c->data), totalSize);
        printf("realloc done \n");
      //  memcpy((void*)(c.data),(void*)(pixels),totalSize);

        for(int i = 0; i < totalSize/3 ; ++i){
          int index = i*3;
          unsigned char B,R,G;
          B = pixels[index];
          G = pixels[index+1];
          R = pixels[index+2];

          c->data[index] = R;
          c->data[index+1] = G;
          c->data[index+2] = B;
        }
        printf("c filled \n");
free(pixels);
   return 1;
}
*/
