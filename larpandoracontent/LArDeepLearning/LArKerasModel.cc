/*Model::KerasModel(const std::string &inputFileName, bool verbose) :
 *
 *  @file   larpandoracontent/LArObjects/LArKerasModel.cc
 *
 *  @brief  Implementation of the lar adaptive boost decision tree class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"
#include "Helpers/XmlHelper.h"

#include "larpandoracontent/LArDeepLearning/LArKerasModel.h"

using namespace pandora;

namespace lar_content
{

KerasModel::KerasModel(const std::string &inputFileName, bool verbose) :
    m_nLayers(std::numeric_limits<int>::max()),
    m_nActiveLayers(std::numeric_limits<int>::max()),
    m_verbose(verbose)
{
    try
    {
        this->LoadWeights(inputFileName);
    }
    catch (const StatusCodeException &statusCodeException)
    {
        std::cout << "KerasModel::LoadWeights - Unable to define network due to unexpected layer type" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::~KerasModel()
{
    for (const Layer *const pLayer : m_layers)
        delete pLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------
/*
StatusCode KerasModel::Initialize(const std::string &bdtXmlFileName, const std::string &bdtName)
{
    return STATUS_CODE_SUCCESS;
}
*/
//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::LoadWeights(const std::string &inputFileName)
{
    if (m_verbose)
        std::cout << "Reading model from " << inputFileName << std::endl;

    std::ifstream inputFileStream(inputFileName.c_str());

    std::string layerType("");
    std::string valueStr("");
    int valueInt(std::numeric_limits<int>::max());

    inputFileStream >> valueStr >> m_nLayers;

    m_nActiveLayers = m_nLayers;

    if (m_verbose)
        std::cout << "nLayers " << m_nLayers << std::endl;

    for (unsigned int layerCount = 0; layerCount < m_nLayers; layerCount++)
    {
        inputFileStream >> valueStr >> valueInt >> layerType;

        if (m_verbose)
            std::cout << "Layer " << valueInt << " " << layerType << ", layerCount = " << layerCount << std::endl;

        Layer *pLayer(nullptr);

        if (layerType == "Conv2D")
        {
            pLayer = new LayerConv2D();
        }
        else if (layerType == "Activation")
        {
           pLayer = new LayerActivation();
        }
        else if (layerType == "MaxPooling2D")
        {
           pLayer = new LayerMaxPooling();
        }
        else if (layerType == "Flatten")
        {
           pLayer = new LayerFlatten();
        }
        else if (layerType == "Dense")
        {
           pLayer = new LayerDense();
        }
        else if (layerType == "Dropout")
        {
           // Dropout layer is not needed in prediction mode
           m_nActiveLayers--;
           continue;
        }

        if (!pLayer)
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pLayer->LoadWeights(inputFileStream);
        m_layers.push_back(pLayer);
    }

    inputFileStream.close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::CalculateOutput(const KerasModel::DataBlock *pDataBlock, Data1D &outputData1D, const pandora::Algorithm *const /*pAlgorithm*/) const
{
    const KerasModel::DataBlock *pInputDataBlock = pDataBlock;
    const KerasModel::DataBlock *pOutputDataBlock(nullptr);

    for (const Layer *const pLayer : m_layers)
    {
        try
        {
            pOutputDataBlock = pLayer->CalculateOutput(pInputDataBlock);
        }
        catch(...)
        {
            std::cout << "KerasModel::CalculateOutput -  Unable to calculate the output for " << pLayer->GetName() << " layer of Keras model" << std::endl;

            // Retain the original data block
            if (pInputDataBlock != pDataBlock)
                delete pInputDataBlock;

            delete pOutputDataBlock;
            return;
        }

        // Retain the original data block
        if (pInputDataBlock != pDataBlock)
            delete pInputDataBlock;

        pInputDataBlock = pOutputDataBlock;
    }

    outputData1D = pOutputDataBlock->GetData1D();
    delete pOutputDataBlock;
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::GetNInputRows() const
{
    if (!m_layers.empty())
        return m_layers.front()->GetNInputRows();

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::GetNInputCols() const
{
    if (!m_layers.empty())
        return m_layers.front()->GetNInputCols();

    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::GetOutputLength() const
{
    int layerCounter(m_layers.size() - 1);
    while (layerCounter > 0 && (m_layers.at(layerCounter)->GetOutputUnits() == 0))
        layerCounter--;

    return m_layers.at(layerCounter)->GetOutputUnits();
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data1D KerasModel::Read1DArray(std::ifstream &inputFileStream, const int nCols)
{
    FloatVector data;
    float value;
    char valueChar;

    inputFileStream >> valueChar; // For [

    for (unsigned int col = 0; col < nCols; col++)
    {
        inputFileStream >> value;
        data.push_back(value);
    }

    inputFileStream >> valueChar; // For ]
    return Data1D(data);
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock::DataBlock()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock::~DataBlock()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::DataBlock::GetDataDim() const
{
    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data1D const &KerasModel::DataBlock::GetData1D() const
{
    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data3D const &KerasModel::DataBlock::GetData3D() const
{
    throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlock::SetData(Data1D const &/*data1D*/)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlock::SetData(Data3D const &/*data3D*/)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlock::ReadFromFile(const std::string &/*inputFileName*/)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock2D::DataBlock2D() :
    m_depth(std::numeric_limits<int>::max()),
    m_rows(std::numeric_limits<int>::max()),
    m_cols(std::numeric_limits<int>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock2D::~DataBlock2D()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::DataBlock2D::GetDataDim() const
{
    return 3;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data3D const &KerasModel::DataBlock2D::GetData3D() const
{
    return m_data3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlock2D::SetData(Data3D const &data3D)
{
    m_data3D = data3D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlock2D::ReadFromFile(const std::string &inputFileName)
{
    std::ifstream inputFile(inputFileName.c_str());
    inputFile >> m_depth >> m_rows >> m_cols;

    for (unsigned int depth = 0; depth < m_depth; depth++)
    {
        Data2D data2D;
        for (unsigned int row = 0; row < m_rows; row++)
        {
            Data1D data1D = KerasModel::Read1DArray(inputFile, m_cols);
            data2D.Append(data1D);
        }
        m_data3D.Append(data2D);
    }
    inputFile.close();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlock2D::ShowName() const
{
    std::cout << "KerasModel::DataBlock2D::ShowName - DataBlock2D " << m_data3D.GetSizeK() << "x" << m_data3D.GetSizeJ() << "x" << m_data3D.GetSizeI() << std::endl;
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlock2D::ShowValues() const
{
    for (unsigned int k = 0; k < m_data3D.GetSizeK(); k++)
    {
        for (unsigned int j = 0; j < m_data3D.GetSizeJ(); j++)
        {
            for (unsigned int i = 0; i < m_data3D.GetSizeI(); i++)
            {
                std::cout << "(" << i << "," << j << "," << k << ") : " << m_data3D.Get(i, j, k) << std::endl;
            }
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlockFlat::DataBlockFlat()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlockFlat::DataBlockFlat(unsigned int size) :
    m_data1D(FloatVector(size))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlockFlat::DataBlockFlat(unsigned int size, float init) :
    m_data1D(FloatVector(size, init))
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlockFlat::~DataBlockFlat()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::DataBlockFlat::GetDataDim() const
{
    return 1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data1D const &KerasModel::DataBlockFlat::GetData1D() const
{
    return m_data1D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlockFlat::SetData(const Data1D &data1D)
{
    m_data1D = data1D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlockFlat::ReadFromFile(const std::string &/*inputFileName*/)
{
    throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlockFlat::ShowName() const
{
    std::cout << "KerasModel::DataBlockFlat::ShowName - DataBlockFlat size " << m_data1D.GetSizeI() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlockFlat::ShowValues() const
{
    std::cout << "KerasModel::DataBlockFlat::ShowValues - DataBlockFlat values : " << std::endl;

    for (int idxI = 0; idxI < m_data1D.GetSizeI(); idxI++)
        std::cout << m_data1D.Get(idxI) << " ";

    std::cout << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::Layer::Layer(const std::string name) :
    m_name(name)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::Layer::~Layer()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

std::string KerasModel::Layer::GetName() const
{
    return m_name;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::LayerFlatten::LayerFlatten() : Layer("Flatten")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::LayerFlatten::LoadWeights(std::ifstream &/*inputFileStream*/)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock* KerasModel::LayerFlatten::CalculateOutput(const KerasModel::DataBlock* pDataBlock) const
{
    const Data3D data3D = pDataBlock->GetData3D();
    KerasModel::DataBlockFlat *pDataBlockFlat = new KerasModel::DataBlockFlat(data3D.GetSizeI() * data3D.GetSizeJ() * data3D.GetSizeK());
    Data1D data1D;

    for (unsigned int k = 0; k < data3D.GetSizeK(); k++)
    {
        for (unsigned int j = 0; j < data3D.GetSizeJ(); j++)
        {
            for (unsigned int i = 0; i < data3D.GetSizeI(); i++)
            {
                data1D.Append(data3D.Get(i,j,k));
            }
        }
    }

    pDataBlockFlat->SetData(data1D);
    return pDataBlockFlat;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerFlatten::GetNInputRows() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerFlatten::GetNInputCols() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerFlatten::GetOutputUnits() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::LayerMaxPooling::LayerMaxPooling() : Layer("MaxPooling2D"),
     m_poolX(std::numeric_limits<int>::max()),
     m_poolY(std::numeric_limits<int>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::LayerMaxPooling::LoadWeights(std::ifstream &inputFileStream)
{
    inputFileStream >> m_poolX >> m_poolY;
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock* KerasModel::LayerMaxPooling::CalculateOutput(const KerasModel::DataBlock* pDataBlock) const
{
    const Data3D data3D = pDataBlock->GetData3D();
    Data3D activeData3D;

    for (unsigned int k = 0; k < data3D.GetSizeK(); k++)
    {
        Data2D activeData2D;
        const unsigned int poolMaxJ(m_poolX != 0 ? data3D.GetSizeJ()/m_poolX : 0);
        const unsigned int poolMaxI(m_poolY != 0 ? data3D.GetSizeI()/m_poolY : 0);

        for (unsigned int j = 0; j < poolMaxJ; j++)
            activeData2D.Append(Data1D(FloatVector(poolMaxI, 0.f)));

        activeData3D.Append(activeData2D);
    }

    for (unsigned int k = 0; k < activeData3D.GetSizeK(); k++)
    {
        for (unsigned int j = 0; j < activeData3D.GetSizeJ(); j++)
        {
            unsigned int startY(j * m_poolY);
            unsigned int endY(startY + m_poolY);

            for (unsigned int i = 0; i < activeData3D.GetSizeI(); i++)
            {
                unsigned int startX(i * m_poolX);
                unsigned int endX(startX + m_poolX);

                FloatVector dataToPool;
                for (unsigned int y = startY; y < endY; y++)
                {
                    for (unsigned int x = startX; x < endX; x++)
                    {
                        dataToPool.push_back(data3D.Get(x, y, k));
                    }
                }

                activeData3D.Set(i, j, k, *std::max_element(dataToPool.begin(), dataToPool.end()));
            }
        }
    }

    KerasModel::DataBlock2D *pDataBlock2D = new KerasModel::DataBlock2D();
    pDataBlock2D->SetData(activeData3D);
    return pDataBlock2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerMaxPooling::GetNInputRows() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerMaxPooling::GetNInputCols() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerMaxPooling::GetOutputUnits() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::LayerActivation::LayerActivation() : Layer("Activation")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::LayerActivation::LoadWeights(std::ifstream &inputFileStream)
{
    inputFileStream >> m_activationType;
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock* KerasModel::LayerActivation::CalculateOutput(const KerasModel::DataBlock* pDataBlock) const
{
    if (pDataBlock->GetDataDim() == 3)
    {
        Data3D data3D = pDataBlock->GetData3D();

        if ("relu" == m_activationType)
        {
            for (unsigned int k = 0; k < data3D.GetSizeK(); k++)
            {
                for (unsigned int j = 0; j < data3D.GetSizeJ(); j++)
                {
                    for (unsigned int i = 0; i < data3D.GetSizeI(); i++)
                    {
                        if(data3D.Get(i, j, k) < 0.f)
                            data3D.Set(i, j, k, 0.f);
                    }
                }
            }

            KerasModel::DataBlock2D *pDataBlock2D = new KerasModel::DataBlock2D();
            pDataBlock2D->SetData(data3D);
            return pDataBlock2D;
        }
        else
        {
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }
    }
    else if (pDataBlock->GetDataDim() == 1)
    {
        Data1D data1D = pDataBlock->GetData1D();

        if ("relu" == m_activationType)
        {
            for (unsigned int i = 0; i < data1D.GetSizeI(); i++)
            {
                if (data1D.Get(i) < 0.f)
                    data1D.Set(i, 0.f);
            }
        }
        else if ("softmax" == m_activationType)
        {
            float sum(0.f);
            for (unsigned int i = 0; i < data1D.GetSizeI(); i++)
            {
                data1D.Set(i, std::exp(data1D.Get(i)));
                sum += data1D.Get(i);
            }
            for (unsigned int i = 0; i < data1D.GetSizeI(); i++)
            {
                data1D.Set(i, data1D.Get(i)/sum);
            }
        }
        else if ("sigmoid" == m_activationType)
        {
            for (unsigned int i = 0; i < data1D.GetSizeI(); i++)
            {
                data1D.Set(i, 1.f / (1.f + std::exp(-1.f * data1D.Get(i))));
            }
        }
        else if ("tanh" == m_activationType)
        {
            for (unsigned int i = 0; i < data1D.GetSizeI(); i++)
            {
                data1D.Set(i, std::tanh(data1D.Get(i)));
            }
        }
        else
        {
            throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }

        KerasModel::DataBlockFlat *pDataBlockFlat = new KerasModel::DataBlockFlat();
        pDataBlockFlat->SetData(data1D);
        return pDataBlockFlat;
    }
    else
    {
        throw StatusCodeException(STATUS_CODE_NOT_ALLOWED);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerActivation::GetNInputRows() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerActivation::GetNInputCols() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerActivation::GetOutputUnits() const
{
    // ATTN: Look at preceeding layer
    return 0;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::LayerConv2D::LayerConv2D() : Layer("Conv2D"),
    m_kernelsCount(std::numeric_limits<int>::max()),
    m_nDepth(std::numeric_limits<int>::max()),
    m_nRows(std::numeric_limits<int>::max()),
    m_nCols(std::numeric_limits<int>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::LayerConv2D::LoadWeights(std::ifstream &inputFileStream)
{
    char valueChar(' ');
    std::string valueStr("");
    float value(std::numeric_limits<float>::max());
    bool skip(false);
    inputFileStream >> m_nRows >> m_nCols >> m_nDepth >> m_kernelsCount >> m_borderMode;

    if (m_borderMode == "[")
    {
        m_borderMode = "valid";
        skip = true;
    }

    for (unsigned int kernel = 0; kernel < m_kernelsCount; kernel++)
    {
        Data3D data3D;
        for (unsigned int depth = 0; depth < m_nDepth; depth++)
        {
            Data2D data2D;
            for (unsigned int row = 0; row < m_nRows; row++)
            {
                if (!skip)
                {
                    // ATTN: For the '['
                    inputFileStream >> valueChar;
                }
                else
                {
                    skip = false;
                }

                Data1D data1D;
                for (unsigned int col = 0; col < m_nCols; col++)
                {
                    inputFileStream >> value;
                    data1D.Append(value);
                }

                // ATTN: For the ']'
                inputFileStream >> valueChar;
                data2D.Append(data1D);
            }
            data3D.Append(data2D);
        }
        m_kernels.Append(data3D);
    }

    // ATTN: For the '['
    inputFileStream >> valueChar;

    for (unsigned int kernel = 0; kernel < m_kernelsCount; kernel++)
    {
        inputFileStream >> value;
        m_bias.Append(value);
    }

    // ATTN: For the '['
    inputFileStream >> valueChar;
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock* KerasModel::LayerConv2D::CalculateOutput(const KerasModel::DataBlock* pDataBlock) const
{
// FIX
    unsigned int startX(std::floor((m_kernels.GetSizeJ() - 1) * 0.5f));
    unsigned int startY(std::floor((m_kernels.GetSizeI() - 1) * 0.5f));

    const Data3D data3D(pDataBlock->GetData3D());
    Data3D activeData3D;

    // Border mode asks whether to shrink grid when concolving or leave it as same size as input
    unsigned int nOutputRows((m_borderMode == "valid") ? data3D.GetSizeJ() - 2 * startX : data3D.GetSizeJ());
    unsigned int nOutputCols((m_borderMode == "valid") ? data3D.GetSizeI() - 2 * startY : data3D.GetSizeI());

    for (unsigned int kernel = 0; kernel < m_kernels.GetSizeL(); kernel++)
    {
        Data2D data2D;

        for(unsigned int rowOutput = 0; rowOutput < nOutputRows; rowOutput++)
        {
            data2D.Append(Data1D(FloatVector(nOutputCols, 0.f)));
        }
        activeData3D.Append(data2D);
    }

    for (unsigned int kernel = 0; kernel < m_kernels.GetSizeL(); kernel++)
    {
        for(unsigned int k = 0; k < data3D.GetSizeK(); k++)
        {
            Data2D data2D(m_borderMode == "valid" ? this->Convolve_SingleDepthValid(data3D.GetData2D(k), m_kernels.GetData2D(k, kernel)) : this->Convolve_SingleDepthSame(data3D.GetData2D(k), m_kernels.GetData2D(k, kernel)));

            for(unsigned int j = 0; j < data2D.GetSizeJ(); j++)
            {
                for(unsigned int i = 0; i < data2D.GetSizeI(); i++)
                {
                    activeData3D.Set(i, j, kernel, activeData3D.Get(i, j, kernel) + data2D.Get(i, j));
                }
            }
        }

        for(unsigned int j = 0; j < activeData3D.GetSizeJ(); j++)
        {
            for(unsigned int i = 0; i < activeData3D.GetSizeI(); i++)
            {
                activeData3D.Set(i, j, kernel, activeData3D.Get(i, j, kernel) + m_bias.Get(kernel));
            }
        }
    }

    KerasModel::DataBlock2D *pDataBlock2D = new KerasModel::DataBlock2D();
    pDataBlock2D->SetData(activeData3D);
    return pDataBlock2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data2D KerasModel::LayerConv2D::Convolve_SingleDepthValid(const Data2D &data2D, const Data2D &filter) const
{
// FIX
    // Function to apply convolution and give output that is reduced in size (e.g. 5x5 * 3x3 = 3x3)
    unsigned int filterSizeJ(filter.GetSizeJ());
    unsigned int filterSizeI(filter.GetSizeI());

    unsigned int startJ((filterSizeJ - 1) * 0.5f);
    unsigned int startI((filterSizeI - 1) * 0.5f);

    // Output data product with correct size to recieve result of convolution
    Data2D activeData2D(data2D.GetSizeJ() - 2 * startJ, Data1D(data2D.GetSizeI() - 2 * startI, 0.f));

    // Loop over all points in the original data set where the filer can be applied
    for (unsigned int j = startJ; j < data2D.GetSizeJ() - startJ; j++)
    {
        for (unsigned int i = startI; i < data2D.GetSizeI() - startI; i++)
        {
            // Apply the filter to the data, seems like the Kernel is flipped so instead of f11*d11 + f12*d12 + ... you'd start with fnn*d11 + fnn-1*d12 + ...
            float sum(0.f);
            for (unsigned int filterJ = 0; filterJ < filterSizeJ; filterJ++)
            {
                for (unsigned int filterI = 0; filterI < filterSizeI; filterI++)
                {
                    sum += filter.Get(filterI, filterJ) * data2D.Get(i - startI + filterI, j - startJ + filterJ);
                }
            }
            activeData2D.Set(i - startI, j - startJ, sum);
        }
    }
    return activeData2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data2D KerasModel::LayerConv2D::Convolve_SingleDepthSame(const Data2D &data2D, const Data2D &filter) const
{
// FIX
    // Function to apply convolution and give output that is the same size as input data
    unsigned int filterSizeJ(filter.GetSizeJ());
    unsigned int filterSizeI(filter.GetSizeI());

    unsigned int startJ((filterSizeJ - 1) * 0.5f);
    unsigned int startI((filterSizeI - 1) * 0.5f);

    unsigned int maxDataJ(data2D.GetSizeJ() - 1);
    unsigned int maxDataI(data2D.GetSizeI() - 1);

    // Output data product with correct size to recieve result of convolution
    Data2D activeData2D(data2D.GetSizeJ(), Data1D(data2D.GetSizeI(), 0.f));

    // Loop over all points in the original data set where the filer can be applied
    for (unsigned int j = startJ; j < data2D.GetSizeJ() - startJ; j++)
    {
        for (unsigned int i = startI; i < data2D.GetSizeI() - startI; i++)
        {
            // Apply the filter to the data, seems like the Kernel is flipped so instead of f11*d11 + f12*d12 + ... you'd start with fnn*d11 + fnn-1*d12 + ...
            float sum(0.f);
            for (unsigned int filterJ = 0; filterJ < filterSizeJ; filterJ++)
            {
                for (unsigned int filterI = 0; filterI < filterSizeI; filterI++)
                {
                    if ((static_cast<int>(j - startJ + filterJ) < 0) || (static_cast<int>(j - startJ + filterJ) > maxDataJ) ||
                        (static_cast<int>(i - startI + filterI) < 0) || (static_cast<int>(i - startI + filterI) > maxDataI))
                        continue;

                    sum += filter.Get(filterI, filterJ) * data2D.Get(i - startI + filterI, j - startJ + filterJ);
                }
            }
            activeData2D.Set(i, j, sum);
        }
    }
    return activeData2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerConv2D::GetNInputRows() const
{
    return m_nRows;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerConv2D::GetNInputCols() const
{
    return m_nCols;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerConv2D::GetOutputUnits() const
{
    return m_kernelsCount;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::LayerDense::LayerDense() : Layer("Dense"),
    m_inputCount(std::numeric_limits<int>::max()),
    m_nNeurons(std::numeric_limits<int>::max())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::LayerDense::LoadWeights(std::ifstream &inputFileStream)
{
    inputFileStream >> m_inputCount >> m_nNeurons;

    char valueChar(' ');
    float value(std::numeric_limits<float>::max());

    for (unsigned int count = 0; count < m_inputCount; count++)
    {
        Data1D data1D;
        inputFileStream >> valueChar; // for '['

        for (unsigned int neuronCount = 0; neuronCount < m_nNeurons; neuronCount++)
        {
            inputFileStream >> value;
            data1D.Append(value);
        }

        inputFileStream >> valueChar; // for ']'
        m_weights.Append(data1D);
    }

    inputFileStream >> valueChar; // for '['

    for (unsigned int neuronCount = 0; neuronCount < m_nNeurons; neuronCount++)
    {
        inputFileStream >> value;
        m_bias.Append(value);
    }

    inputFileStream >> valueChar; // for ']'
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock* KerasModel::LayerDense::CalculateOutput(const KerasModel::DataBlock* pDataBlock) const
{
    KerasModel::DataBlockFlat *pDataBlockFlat = new KerasModel::DataBlockFlat(m_nNeurons, 0.f);
    Data1D activeData1D(pDataBlockFlat->GetData1D());
    const Data1D data1D(pDataBlock->GetData1D());

    for (unsigned int weightCount = 0; weightCount < m_weights.GetSizeJ(); weightCount++)
    {
        const Data1D weights(m_weights.GetData1D(weightCount));
        float data(data1D.Get(weightCount));

        unsigned int neuronCounter(0);

        while (neuronCounter < m_nNeurons)
        {
            activeData1D.Set(neuronCounter, activeData1D.Get(neuronCounter) + (weights.Get(neuronCounter) * data));
            neuronCounter++;
        }
    }

    for (unsigned int biasCount = 0; biasCount < m_nNeurons; biasCount++)
    {
        activeData1D.Set(biasCount, activeData1D.Get(biasCount) + m_bias.Get(biasCount));
    }

    pDataBlockFlat->SetData(activeData1D);
    return pDataBlockFlat;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerDense::GetNInputRows() const
{
    return 1;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerDense::GetNInputCols() const
{
    return m_inputCount;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerDense::GetOutputUnits() const
{
    return m_nNeurons;
}

} // namespace lar_content
