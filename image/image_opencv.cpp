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

    OpenGLMatrix modelv;
    convertOpenCVProjectionMatrixToOpenGL(c->cameraMatrix,0.1,3,c->projM.m);
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
