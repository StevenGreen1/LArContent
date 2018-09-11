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

KerasModel::KerasModel()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::~KerasModel()
{
    for (const Layer *const pLayer : m_layers)
        delete pLayer;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KerasModel::Initialize(const std::string &cnnXmlFileName, const std::string &cnnName)
{
    TiXmlDocument xmlDocument(cnnXmlFileName);

    if (!xmlDocument.LoadFile())
    {
        std::cout << "KerasModel::Initialize - Invalid xml file." << std::endl;
        return STATUS_CODE_INVALID_PARAMETER;
    }

    const TiXmlHandle xmlDocumentHandle(&xmlDocument);
    TiXmlNode *pContainerXmlNode(TiXmlHandle(xmlDocumentHandle).FirstChildElement().Element());

    while (pContainerXmlNode)
    {
        if (pContainerXmlNode->ValueStr() != "ConvolutionalNeuralNetwork")
            return STATUS_CODE_FAILURE;

        const TiXmlHandle currentHandle(pContainerXmlNode);

        std::string currentName;
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(currentHandle, "Name", currentName));

        if (currentName.empty() || (currentName.size() > 1000))
        {
            std::cout << "KerasModel::Initialize - Implausible ConvolutionalNeuralNetwork name extracted from xml." << std::endl;
            return STATUS_CODE_INVALID_PARAMETER;
        }

        if (currentName == cnnName)
            break;

        pContainerXmlNode = pContainerXmlNode->NextSibling();
    }

    if (!pContainerXmlNode)
    {
        std::cout << "AdaBoostDecisionTree: Could not find an ConvolutionalNeuralNetwork of name " << cnnName << std::endl;
        return STATUS_CODE_NOT_FOUND;
    }

    const TiXmlHandle xmlHandle(pContainerXmlNode);
    this->LoadWeights(&xmlHandle);
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::LoadWeights(const TiXmlHandle *const pXmlHandle)
{
    TiXmlElement *pCurrentXmlElement = pXmlHandle->FirstChild().Element();

    while (pCurrentXmlElement)
    {
        if (STATUS_CODE_SUCCESS != this->ReadComponent(pCurrentXmlElement))
            throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);

        pCurrentXmlElement = pCurrentXmlElement->NextSiblingElement();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode KerasModel::ReadComponent(TiXmlElement *pCurrentXmlElement)
{
    const std::string componentName(pCurrentXmlElement->ValueStr());
    TiXmlHandle currentHandle(pCurrentXmlElement);

    if ((std::string("Name") == componentName) || (std::string("Timestamp") == componentName))
        return STATUS_CODE_SUCCESS;

    if (std::string("LayerConv2D") == componentName)
    {
        m_layers.emplace_back(new LayerConv2D(&currentHandle));
        return STATUS_CODE_SUCCESS;
    }

    if (std::string("LayerActivation") == componentName)
    {
        m_layers.emplace_back(new LayerActivation(&currentHandle));
        return STATUS_CODE_SUCCESS;
    }

    if (std::string("LayerMaxPooling2D") == componentName)
    {
        m_layers.emplace_back(new LayerMaxPooling(&currentHandle));
        return STATUS_CODE_SUCCESS;
    }

    if (std::string("LayerFlatten") == componentName)
    {
        m_layers.emplace_back(new LayerFlatten(&currentHandle));
        return STATUS_CODE_SUCCESS;
    }

    if (std::string("LayerDense") == componentName)
    {
        m_layers.emplace_back(new LayerDense(&currentHandle));
        return STATUS_CODE_SUCCESS;
    }

    return STATUS_CODE_INVALID_PARAMETER;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::CalculateOutput(const KerasModel::DataBlock *pDataBlock, Data1D &outputData1D, const pandora::Algorithm *const pAlgorithm) const
{
    const KerasModel::DataBlock *pInputDataBlock = pDataBlock;
    const KerasModel::DataBlock *pOutputDataBlock(nullptr);

    for (const Layer *const pLayer : m_layers)
    {
        try
        {
            pOutputDataBlock = pLayer->CalculateOutput(pInputDataBlock);

            if (pLayer->GetName() == "Conv2D")
            {
                const int m_gridSize(14);
                Data3D outputData3D = pOutputDataBlock->GetData3D();
                for (unsigned int k = 0; k < outputData3D.GetSizeK(); k++)
                {
                    TwoDHistogram twoDHistogram(m_gridSize, 0.f,  m_gridSize, m_gridSize, 0.f,  m_gridSize);
                    for (unsigned int j = 0; j < outputData3D.GetSizeJ(); j++)
                    {
                        for (unsigned int i = 0; i < outputData3D.GetSizeI(); i++)
                        {
                            twoDHistogram.Fill(i,j,outputData3D.Get(i,j,k));
                        }
                    }
                    PandoraMonitoringApi::DrawPandoraHistogram(pAlgorithm->GetPandora(), twoDHistogram, "COLZ");
                }
            }
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
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock2D::DataBlock2D()
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

void KerasModel::DataBlockFlat::ShowName() const
{
    std::cout << "KerasModel::DataBlockFlat::ShowName - DataBlockFlat size " << m_data1D.GetSizeI() << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void KerasModel::DataBlockFlat::ShowValues() const
{
    std::cout << "KerasModel::DataBlockFlat::ShowValues - DataBlockFlat values : ";

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

KerasModel::LayerFlatten::LayerFlatten(const TiXmlHandle *const /*pXmlHandle*/) : Layer("Flatten")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock* KerasModel::LayerFlatten::CalculateOutput(const KerasModel::DataBlock* pDataBlock) const
{
    const Data3D data3D = pDataBlock->GetData3D();
    KerasModel::DataBlockFlat *pDataBlockFlat = new KerasModel::DataBlockFlat(data3D.GetSizeI() * data3D.GetSizeJ() * data3D.GetSizeK());
    Data1D data1D;

    // ATTN: This ordering is expected by downstream layers, do not change
    for (unsigned int j = 0; j < data3D.GetSizeJ(); j++)
    {
        for (unsigned int i = 0; i < data3D.GetSizeI(); i++)
        {
            for (unsigned int k = 0; k < data3D.GetSizeK(); k++)
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

KerasModel::LayerMaxPooling::LayerMaxPooling(const TiXmlHandle *const pXmlHandle) : Layer("MaxPooling2D"),
     m_poolX(std::numeric_limits<int>::max()),
     m_poolY(std::numeric_limits<int>::max())
{
    for (TiXmlElement *pHeadTiXmlElement = pXmlHandle->FirstChildElement().ToElement(); pHeadTiXmlElement != NULL; pHeadTiXmlElement = pHeadTiXmlElement->NextSiblingElement())
    {
        if ("MaxPooling2DConfig" == pHeadTiXmlElement->ValueStr())
        {
            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("PoolX", &m_poolX))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("PoolY", &m_poolY))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }
    }
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
        {
            activeData2D.Append(Data1D(FloatVector(poolMaxI, 0.f)));
        }

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

KerasModel::LayerActivation::LayerActivation(const TiXmlHandle *const pXmlHandle) : Layer("Activation")
{
    pandora::StringVector activationTypes = {"relu", "softmax", "tanh", "sigmoid"};

    for (TiXmlElement *pHeadTiXmlElement = pXmlHandle->FirstChildElement().ToElement(); pHeadTiXmlElement != NULL; pHeadTiXmlElement = pHeadTiXmlElement->NextSiblingElement())
    {
        if ("ActivationFunction" == pHeadTiXmlElement->ValueStr())
        {
            m_activationType = pHeadTiXmlElement->GetText();

            if (std::find(activationTypes.begin(), activationTypes.end(), m_activationType) == activationTypes.end())
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);
        }
    }
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

KerasModel::LayerConv2D::LayerConv2D(const TiXmlHandle *const pXmlHandle) : Layer("Conv2D"),
    m_nKernels(std::numeric_limits<int>::max()),
    m_nDeep(std::numeric_limits<int>::max()),
    m_nRows(std::numeric_limits<int>::max()),
    m_nCols(std::numeric_limits<int>::max())
{
    bool configLoaded(false);

    for (TiXmlElement *pHeadTiXmlElement = pXmlHandle->FirstChildElement().ToElement(); pHeadTiXmlElement != NULL; pHeadTiXmlElement = pHeadTiXmlElement->NextSiblingElement())
    {
        if ("LayerConv2DConfig" == pHeadTiXmlElement->ValueStr())
        {
            m_borderMode = pHeadTiXmlElement->Attribute("borderMode");

            if (m_borderMode.empty())
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("nKernels", &m_nKernels))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("nDeep", &m_nDeep))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("nRows", &m_nRows))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("nCols", &m_nCols))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            m_kernels = Data4D(m_nKernels, Data3D(m_nDeep, Data2D(m_nRows, Data1D(m_nCols, std::numeric_limits<float>::max()))));
            m_bias = Data1D(m_nKernels, 0.f);
            configLoaded = true;
        }
        else if ("Conv2DWeight" == pHeadTiXmlElement->ValueStr())
        {
            if (!configLoaded)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            int kernel(std::numeric_limits<int>::max()), depth(std::numeric_limits<int>::max()), row(std::numeric_limits<int>::max()), col(std::numeric_limits<int>::max());
            double value(std::numeric_limits<double>::max());

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("kernel", &kernel))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("depth", &depth))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("row", &row))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("col", &col))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryDoubleAttribute("value", &value))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            m_kernels.Set(col, row, depth, kernel, value);
        }
        else if ("Conv2DBias" == pHeadTiXmlElement->ValueStr())
        {
            if (!configLoaded)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            int kernel(std::numeric_limits<int>::max());
            double value(std::numeric_limits<double>::max());

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("kernel", &kernel))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryDoubleAttribute("bias", &value))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            m_bias.Set(kernel, value);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock* KerasModel::LayerConv2D::CalculateOutput(const KerasModel::DataBlock* pDataBlock) const
{
    unsigned int startX(std::floor((m_kernels.GetSizeJ() - 1) * 0.5f));
    unsigned int startY(std::floor((m_kernels.GetSizeI() - 1) * 0.5f));

    const Data3D data3D(pDataBlock->GetData3D());
    Data3D activeData3D;

    unsigned int nOutputRows((m_borderMode == "valid") ? data3D.GetSizeJ() - 2 * startX : data3D.GetSizeJ());
    unsigned int nOutputCols((m_borderMode == "valid") ? data3D.GetSizeI() - 2 * startY : data3D.GetSizeI());

    for (unsigned int l = 0; l < m_kernels.GetSizeL(); l++)
    {
        Data2D data2D;

        for(unsigned int rowOutput = 0; rowOutput < nOutputRows; rowOutput++)
        {
            data2D.Append(Data1D(FloatVector(nOutputCols, 0.f)));
        }
        activeData3D.Append(data2D);
    }

    for (unsigned int l = 0; l < m_kernels.GetSizeL(); l++)
    {
        for(unsigned int k = 0; k < data3D.GetSizeK(); k++)
        {
            Data2D data2D(m_borderMode == "valid" ? this->Convolve_SingleDepthValid(data3D.GetData2D(k), m_kernels.GetData2D(k, l)) : this->Convolve_SingleDepthSame(data3D.GetData2D(k), m_kernels.GetData2D(k, l)));

            for(unsigned int j = 0; j < data2D.GetSizeJ(); j++)
            {
                for(unsigned int i = 0; i < data2D.GetSizeI(); i++)
                {
                    activeData3D.Set(i, j, l, activeData3D.Get(i, j, l) + data2D.Get(i, j));
                }
            }
        }

        for(unsigned int j = 0; j < activeData3D.GetSizeJ(); j++)
        {
            for(unsigned int i = 0; i < activeData3D.GetSizeI(); i++)
            {
                activeData3D.Set(i, j, l, activeData3D.Get(i, j, l) + m_bias.Get(l));
            }
        }
    }

    KerasModel::DataBlock2D *pDataBlock2D = new KerasModel::DataBlock2D();
    pDataBlock2D->SetData(activeData3D);
    pDataBlock2D->ShowValues();
    return pDataBlock2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data2D KerasModel::LayerConv2D::Convolve_SingleDepthValid(const Data2D &data2D, const Data2D &kernel) const
{
    unsigned int kernelSizeJ(kernel.GetSizeJ());
    unsigned int kernelSizeI(kernel.GetSizeI());

    unsigned int startJ((kernelSizeJ - 1) * 0.5f);
    unsigned int startI((kernelSizeI - 1) * 0.5f);

    Data2D activeData2D(data2D.GetSizeJ() - 2 * startJ, Data1D(data2D.GetSizeI() - 2 * startI, 0.f));

    for (unsigned int j = startJ; j < data2D.GetSizeJ() - startJ; j++)
    {
        for (unsigned int i = startI; i < data2D.GetSizeI() - startI; i++)
        {
            float sum(0.f);
            for (unsigned int kernelJ = 0; kernelJ < kernelSizeJ; kernelJ++)
            {
                for (unsigned int kernelI = 0; kernelI < kernelSizeI; kernelI++)
                {
                    sum += kernel.Get(kernelI, kernelJ) * data2D.Get(i - startI + kernelI, j - startJ + kernelJ);
                }
            }
            activeData2D.Set(i - startI, j - startJ, sum);
        }
    }
    return activeData2D;
}

//------------------------------------------------------------------------------------------------------------------------------------------

Data2D KerasModel::LayerConv2D::Convolve_SingleDepthSame(const Data2D &data2D, const Data2D &kernel) const
{
    unsigned int kernelSizeJ(kernel.GetSizeJ());
    unsigned int kernelSizeI(kernel.GetSizeI());

    unsigned int startJ((kernelSizeJ - 1) * 0.5f);
    unsigned int startI((kernelSizeI - 1) * 0.5f);

    unsigned int maxDataJ(data2D.GetSizeJ() - 1);
    unsigned int maxDataI(data2D.GetSizeI() - 1);

    Data2D activeData2D(data2D.GetSizeJ(), Data1D(data2D.GetSizeI(), 0.f));

    for (unsigned int j = startJ; j < data2D.GetSizeJ() - startJ; j++)
    {
        for (unsigned int i = startI; i < data2D.GetSizeI() - startI; i++)
        {
            float sum(0.f);
            for (unsigned int kernelJ = 0; kernelJ < kernelSizeJ; kernelJ++)
            {
                for (unsigned int kernelI = 0; kernelI < kernelSizeI; kernelI++)
                {
                    if ((static_cast<int>(j - startJ + kernelJ) < 0) || (static_cast<int>(j - startJ + kernelJ) > maxDataJ) ||
                        (static_cast<int>(i - startI + kernelI) < 0) || (static_cast<int>(i - startI + kernelI) > maxDataI))
                        continue;

                    sum += kernel.Get(kernelI, kernelJ) * data2D.Get(i - startI + kernelI, j - startJ + kernelJ);
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
    return m_nKernels;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::LayerDense::LayerDense(const TiXmlHandle *const pXmlHandle) : Layer("Dense"),
    m_nInputNodes(std::numeric_limits<int>::max()),
    m_nOutputNodes(std::numeric_limits<int>::max())
{
    bool configLoaded(false);

    for (TiXmlElement *pHeadTiXmlElement = pXmlHandle->FirstChildElement().ToElement(); pHeadTiXmlElement != NULL; pHeadTiXmlElement = pHeadTiXmlElement->NextSiblingElement())
    {
        if ("LayerDenseConfig" == pHeadTiXmlElement->ValueStr())
        {
            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("nInputNodes", &m_nInputNodes))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("nOutputNodes", &m_nOutputNodes))
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            m_weights = Data2D(m_nInputNodes, Data1D(m_nOutputNodes, std::numeric_limits<float>::max()));
            m_bias = Data1D(m_nInputNodes, std::numeric_limits<float>::max());
            configLoaded = true;
        }
        else if ("DenseWeight" == pHeadTiXmlElement->ValueStr())
        {
            if (!configLoaded)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            int connection(std::numeric_limits<int>::max());

            if (TIXML_SUCCESS != pHeadTiXmlElement->QueryIntAttribute("connection", &connection))
                    throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            for (int i = 0; i < m_nOutputNodes; i++)
            {
                double value(std::numeric_limits<double>::max());
                std::string name("weight" + std::to_string(i));

                if (TIXML_SUCCESS != pHeadTiXmlElement->QueryDoubleAttribute(name.c_str(), &value))
                    throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                m_weights.Set(i, connection, value);
            }
        }
        else if ("DenseBias" == pHeadTiXmlElement->ValueStr())
        {
            if (!configLoaded)
                throw StatusCodeException(STATUS_CODE_NOT_FOUND);

            for (int i = 0; i < m_nOutputNodes; i++)
            {
                double value(std::numeric_limits<double>::max());
                std::string name("bias" + std::to_string(i));

                if (TIXML_SUCCESS != pHeadTiXmlElement->QueryDoubleAttribute(name.c_str(), &value))
                    throw StatusCodeException(STATUS_CODE_NOT_FOUND);

                m_bias.Set(i, value);
            }
        }
    }

    if (m_weights.GetSizeJ() != m_nInputNodes || m_bias.GetSizeI() != m_nInputNodes)
        throw StatusCodeException(STATUS_CODE_NOT_FOUND);
}

//------------------------------------------------------------------------------------------------------------------------------------------

KerasModel::DataBlock* KerasModel::LayerDense::CalculateOutput(const KerasModel::DataBlock* pDataBlock) const
{
    KerasModel::DataBlockFlat *pDataBlockFlat = new KerasModel::DataBlockFlat(m_nOutputNodes, 0.f);
    Data1D activeData1D(pDataBlockFlat->GetData1D());
    const Data1D data1D(pDataBlock->GetData1D());

    for (unsigned int j = 0; j < m_weights.GetSizeJ(); j++)
    {
        const Data1D connectionWeights(m_weights.GetData1D(j));
        float inputData(data1D.Get(j));

        unsigned int i(0);

        while (i < m_nOutputNodes)
        {
            activeData1D.Set(i, activeData1D.Get(i) + (connectionWeights.Get(i) * inputData));
            i++;
        }
    }

    for (unsigned int i = 0; i < m_nOutputNodes; i++)
    {
        activeData1D.Set(i, activeData1D.Get(i) + m_bias.Get(i));
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
    return m_nInputNodes;
}

//------------------------------------------------------------------------------------------------------------------------------------------

unsigned int KerasModel::LayerDense::GetOutputUnits() const
{
    return m_nOutputNodes;
}

} // namespace lar_content
