#include "IO.h"



void screenshot_bmp(char *fileName, unsigned int width, unsigned int height) {
    const size_t format_nchannels = 3;
    GLubyte *pixels = (GLubyte *) malloc(format_nchannels * sizeof(GLubyte) * width * height);
    glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
    int32 bytesPerPixel = 3;
    //Open file in binary mode
    FILE *outputFile = fopen(fileName, "wb");
    //*****HEADER************//
    //write signature
    const char *BM = "BM";
    fwrite(&BM[0], 1, 1, outputFile);
    fwrite(&BM[1], 1, 1, outputFile);
    //Write file size considering padded bytes
    int paddedRowSize = (int)(4 * ceil((float)width/4.0f))*bytesPerPixel;
    int32 fileSize = paddedRowSize*height + HEADER_SIZE + INFO_HEADER_SIZE;
    fwrite(&fileSize, 4, 1, outputFile);
    //Write reserved
    int32 reserved = 0x0000;
    fwrite(&reserved, 4, 1, outputFile);
    //Write data offset
    int32 dataOffset = HEADER_SIZE+INFO_HEADER_SIZE;
    fwrite(&dataOffset, 4, 1, outputFile);

    //*******INFO*HEADER******//
    //Write size
    int32 infoHeaderSize = INFO_HEADER_SIZE;
    fwrite(&infoHeaderSize, 4, 1, outputFile);
    //Write width and height
    fwrite(&width, 4, 1, outputFile);
    fwrite(&height, 4, 1, outputFile);
    //Write planes
    int16 planes = 1; //always 1
    fwrite(&planes, 2, 1, outputFile);
    //write bits per pixel
    int16 bitsPerPixel = bytesPerPixel * 8;
    fwrite(&bitsPerPixel, 2, 1, outputFile);
    //write compression
    int32 compression = NO_COMPRESION;
    fwrite(&compression, 4, 1, outputFile);
    //write image size (in bytes)
    int32 imageSize = width*height*bytesPerPixel;
    fwrite(&imageSize, 4, 1, outputFile);
    //write resolution (in pixels per meter)
    int32 resolutionX = 11811; //300 dpi
    int32 resolutionY = 11811; //300 dpi
    fwrite(&resolutionX, 4, 1, outputFile);
    fwrite(&resolutionY, 4, 1, outputFile);
    //write colors used
    int32 colorsUsed = MAX_NUMBER_OF_COLORS;
    fwrite(&colorsUsed, 4, 1, outputFile);
    //Write important colors
    int32 importantColors = ALL_COLORS_REQUIRED;
    fwrite(&importantColors, 4, 1, outputFile);

    int unpaddedRowSize = width*bytesPerPixel;
    //conversion form rgb to bgr
    int totalSize = unpaddedRowSize*(height);
    GLubyte *pixels_up = (GLubyte *) malloc(format_nchannels * sizeof(GLubyte) * width * height);
    for(int i = 0; i < totalSize/3 ; ++i){
      int index = i*3;
      unsigned char B,R,G;
      B = pixels[index];
      G = pixels[index+1];
      R = pixels[index+2];

      pixels_up[index] = R;
      pixels_up[index+1] = G;
      pixels_up[index+2] = B;
    }

free(pixels);
    //write data
    int i = 0;

    for ( i = 0; i < height; i++)
    {
            //start writing from the beginning of last row in the pixel array
            int pixelOffset = i*unpaddedRowSize;//((height - i) - 1)*unpaddedRowSize;
            fwrite(&pixels_up[pixelOffset], 1, paddedRowSize, outputFile);
    }
    fclose(outputFile);
    free(pixels_up);
}
