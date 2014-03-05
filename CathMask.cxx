/*=========================================================================

  Program:   3D Slicer Catheter Mask CLI
  Language:  C++
  Author:    Junichi Tokuda, Ph.D. (Brigham and Women's Hospital)

  Copyright (c) Brigham and Women's Hospital. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"

#include "itkTransformFileReader.h"
#include "itkEllipseSpatialObject.h"
#include "itkSpatialObjectToImageFilter.h"
#include "itkGroupSpatialObject.h"
#include "itkAffineTransform.h"

#include "itkPluginUtilities.h"
#include "CathMaskCLP.h"


// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace {

template<class T> int DoIt( int argc, char * argv[], T )
{
  PARSE_ARGS;

  const     unsigned int        Dimension       = 3;

  typedef   T                   FileInputPixelType;
  typedef   float               InternalPixelType;
  typedef   int                 OutputPixelType;

  typedef   itk::Image< FileInputPixelType, Dimension > FileInputImageType;
  typedef   itk::Image< InternalPixelType, Dimension >  InternalImageType;
  typedef   itk::Image< OutputPixelType, Dimension >    OutputImageType;
  typedef   itk::Image< OutputPixelType, Dimension >    LabelImageType;

  typedef   itk::ImageFileReader< InternalImageType >  ReaderType;
  typename ReaderType::Pointer reader = ReaderType::New();  

  // Read reference image (just to get size information)
  reader->SetFileName( referenceVolume.c_str() );
  reader->Update();

  // Read transform
  typedef itk::TransformFileReader TransformReaderType;
  TransformReaderType::TransformType::Pointer cathTransform;
  TransformReaderType::Pointer cathTransformReader;  

  typedef itk::AffineTransform< double, 3 > AffineTransformType;
  AffineTransformType::Pointer atrans = AffineTransformType::New();
  atrans->SetIdentity();

  if( catheterTransform != "" )
    {
    cathTransformReader = TransformReaderType::New();
    cathTransformReader->SetFileName( catheterTransform );
    
    try
      {
      cathTransformReader->Update();
      }
    catch( itk::ExceptionObject & err )
      {
      std::cerr << err << std::endl;
      return EXIT_FAILURE;
      }
    
    if (cathTransformReader->GetTransformList()->size() != 0 )
      {
      cathTransform = *(cathTransformReader->GetTransformList()->begin() );
      }

    atrans->SetParameters(cathTransform->GetParameters());
    std::cerr << "transform: " << atrans->GetMatrix() << std::endl;
    }

  AffineTransformType::InverseTransformBaseType::Pointer itrans;
  itrans = atrans->GetInverseTransform();

  AffineTransformType::Pointer trans;
  trans = dynamic_cast<AffineTransformType *>(itrans.GetPointer() );

  //AffineTransformType::MatrixType matrix;
  //AffineTransformType::OffsetType offset;
  //matrix = trans->GetMatrix();
  //offset = trans->GetOffset();

  // Generate ellipse
  typedef itk::GroupSpatialObject< Dimension >     GroupType;

  typedef itk::EllipseSpatialObject< Dimension > EllipseType;
  //typedef itk::EllipseSpatialObject EllipseType;
  EllipseType::Pointer ellipse = EllipseType::New();

  
  EllipseType::ArrayType dim;
  dim[0] = diameter/2.0;
  dim[1] = diameter/2.0;
  dim[2] = length/2.0;
  ellipse->SetRadius(dim);
  //ellipse->SetRadius(diameter/2.0);
  //ellipse->SetHeight(length);
  
  typedef GroupType::TransformType  GroupTransformType;

  //GroupTransformType::Pointer transform = GroupTransformType::New();
  //transform->SetIdentity();
  //GroupTransformType::OutputVectorType  translation;
  //GroupTransformType::CenterType        center;
  //translation[ 0 ] =  offset[0];
  //translation[ 1 ] =  offset[1];
  //translation[ 2 ] =  offset[2];
  //
  //transform->Translate( translation, false );
  //
  //ellipse->SetObjectToParentTransform( transform );
  typedef itk::ScalableAffineTransform< double, 3 > ScalableAffineTransformType;
  ScalableAffineTransformType::Pointer ellipseTransform = ScalableAffineTransformType::New();
  ellipseTransform->SetIdentity();
  ellipseTransform->SetMatrix(trans->GetMatrix());
  ellipseTransform->SetOffset(trans->GetOffset());
  ellipse->SetObjectToParentTransform( ellipseTransform );
  ellipse->ComputeObjectToWorldTransform();
  //ellipse->ComputeLocalBoundingBox();
  //ellipse->ComputeBoundingBox();

  
  //GroupType::Pointer group = GroupType::New();
  //group->AddSpatialObject( ellipse );
  //std::cerr << group->GetBoundingBox() << std::endl;

  //typedef itk::SpatialObjectToImageFilter< GroupType, LabelImageType >   SpatialObjectToImageFilterType;
  typedef itk::SpatialObjectToImageFilter< EllipseType, LabelImageType >   SpatialObjectToImageFilterType;
  SpatialObjectToImageFilterType::Pointer imageFilter = SpatialObjectToImageFilterType::New();  

  // Copy image dimensions from the reference image
  imageFilter->SetDirection( reader->GetOutput()->GetDirection() );
  imageFilter->SetOrigin( reader->GetOutput()->GetOrigin() );
  imageFilter->SetSpacing( reader->GetOutput()->GetSpacing() );
  imageFilter->SetSize( reader->GetOutput()->GetLargestPossibleRegion().GetSize() );
  //imageFilter->SetInput( group );
  imageFilter->SetInput( ellipse );
  //std::cerr << imageFilter->GetBoundingBox() << std::endl;

  const OutputPixelType insideUnits  = 1;
  const OutputPixelType outsideUnits = 0;

  ellipse->SetDefaultInsideValue( insideUnits );
  ellipse->SetDefaultOutsideValue( outsideUnits );
  imageFilter->SetUseObjectValue( true );
  imageFilter->SetOutsideValue( outsideUnits );
  imageFilter->Update();

  // Output image
  typedef itk::ImageFileWriter< OutputImageType >     WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( outputVolume.c_str() );
  writer->SetInput( imageFilter->GetOutput() );
  writer->SetUseCompression(1);

  try
    {
    writer->Update();
    }
  catch (itk::ExceptionObject &err)
    {
    std::cerr << err << std::endl;
    return EXIT_FAILURE ;
    }

  return EXIT_SUCCESS;
}

} // end of anonymous namespace


int main( int argc, char * argv[] )
{
  PARSE_ARGS;

  itk::ImageIOBase::IOPixelType pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    itk::GetImageType (referenceVolume, pixelType, componentType);

    // This filter handles all types on input, but only produces
    // signed types
    switch (componentType)
      {
      case itk::ImageIOBase::UCHAR:
        return DoIt( argc, argv, static_cast<unsigned char>(0));
        break;
      case itk::ImageIOBase::CHAR:
        return DoIt( argc, argv, static_cast<char>(0));
        break;
      case itk::ImageIOBase::USHORT:
        return DoIt( argc, argv, static_cast<unsigned short>(0));
        break;
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0));
        break;
      case itk::ImageIOBase::UINT:
        return DoIt( argc, argv, static_cast<unsigned int>(0));
        break;
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0));
        break;
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0));
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0));
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0));
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<double>(0));
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }

  catch( itk::ExceptionObject &excep)
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
