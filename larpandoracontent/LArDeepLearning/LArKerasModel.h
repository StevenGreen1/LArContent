/**
 *  @file   larpandoracontent/LArDeepLearning/LArKerasModel.h
 *
 *  @brief  Header file for the lar keras model class.
 *
 *  $Log: $
 */
#ifndef LAR_KERAS_MODEL
#define LAR_KERAS_MODEL 1

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArObjects/LArMvaInterface.h"

#include "Pandora/StatusCodes.h"

#include <functional>
#include <map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  KerasModel class
 */
class KerasModel
{
public:
    typedef pandora::FloatVector Data1D;
    typedef std::vector<Data1D> Data2D;
    typedef std::vector<Data2D> Data3D;
    typedef std::vector<Data3D> Data4D;

    /**
     *  @brief  Constructor.
     */
    KerasModel(const std::string &inputFileName, bool verbose);

    /**
     *  @brief  Destructor.
     */
    ~KerasModel();

    //pandora::StatusCode Initialize(const std::string &parameterLocation, const std::string &bdtName);

    class DataBlock
    {
    public:
        /**
         *  @brief  Constructor.
         */
        DataBlock();

        /**
         *  @brief  Destructor.
         */
        virtual ~DataBlock();

        virtual unsigned int GetDataDim() const;

        // ATTN: Returning const value isn't a sensible thing to do.  Change this.
        virtual Data1D const &GetData1D() const;

        virtual Data3D const &GetData3D() const;

        virtual void SetData(const Data1D &data1D);

        virtual void SetData(const Data3D &data3D);

        virtual void ReadFromFile(const std::string &inputFileName);

        virtual void ShowName() const = 0;

        virtual void ShowValues() const = 0;
    };

    class DataBlock2D : public DataBlock
    {
    public:
        DataBlock2D();

        ~DataBlock2D();

        unsigned int GetDataDim() const;

        Data3D const &GetData3D() const;

        void SetData(const Data3D &data3D);

        void ReadFromFile(const std::string &inputFileName);

        void ShowName() const;

        void ShowValues() const;

    private:
        Data3D     m_data3D;  ///< Data
        int        m_depth;   ///< Depth of data
        int        m_rows;    ///< Number of rows in data
        int        m_cols;    ///< Number of columns in data
    };

    class DataBlockFlat : public DataBlock
    {
    public:
        DataBlockFlat();

        DataBlockFlat(unsigned int size);

        DataBlockFlat(unsigned int size, float init);

        ~DataBlockFlat();

        unsigned int GetDataDim() const;

        Data1D const &GetData1D() const;

        void SetData(const Data1D &data1D);

        void ReadFromFile(const std::string &inputFileName);

        void ShowName() const;

        void ShowValues() const;

    private:
        Data1D     m_data1D;   ///< Data
    };

    class Layer
    {
    public:
        Layer(const std::string name);

        virtual ~Layer();

        virtual void LoadWeights(std::ifstream &inputFileStream) = 0;

        virtual DataBlock* CalculateOutput(DataBlock* pDataBlock) = 0;

        virtual unsigned int GetNInputRows() const = 0;

        virtual unsigned int GetNInputCols() const = 0;

        virtual unsigned int GetOutputUnits() const = 0;

        std::string GetName() const;

    private:
        std::string     m_name;    ///< Layer name
    };

    class LayerFlatten : public Layer
    {
    public:
        LayerFlatten();

        void LoadWeights(std::ifstream &inputFileStream);

        DataBlock* CalculateOutput(DataBlock* pDataBlock);

        unsigned int GetNInputRows() const;

        unsigned int GetNInputCols() const;

        unsigned int GetOutputUnits() const;
    };

    class LayerMaxPooling : public Layer
    {
    public:
        LayerMaxPooling();

        void LoadWeights(std::ifstream &inputFileStream);

        DataBlock* CalculateOutput(DataBlock* pDataBlock);

        unsigned int GetNInputRows() const;

        unsigned int GetNInputCols() const;

        unsigned int GetOutputUnits() const;

    private:
        int     m_poolX;    ///< The pool size along x
        int     m_poolY;    ///< The pool size along y
    };

    class LayerActivation : public Layer
    {
    public:
        LayerActivation();

        void LoadWeights(std::ifstream &inputFileStream);

        DataBlock* CalculateOutput(DataBlock* pDataBlock);

        unsigned int GetNInputRows() const;

        unsigned int GetNInputCols() const;

        unsigned int GetOutputUnits() const;

    private:
        std::string     m_activationType;    ///< Activation type
    };

    class LayerConv2D : public Layer
    {
    public:
        LayerConv2D();

        void LoadWeights(std::ifstream &inputFileStream);

        DataBlock* CalculateOutput(DataBlock* pDataBlock);

        unsigned int GetNInputRows() const;

        unsigned int GetNInputCols() const;

        unsigned int GetOutputUnits() const;

        /**
         *  @brief  Convolve data with a filter and give output that is of reduced size
         *
         *  @param  data2D input image
         *  @param  filter to apply to image
         *
         *  @return The convolved image
         */
        Data2D Convolve_SingleDepthValid(const Data2D &data2D, const Data2D &filter) const;

        /**
         *  @brief  Convolve data with a filter and give output that is of the same size as the input
         *
         *  @param  data2D input image
         *  @param  filter to apply to image
         *
         *  @return The convolved image
         */
        Data2D Convolve_SingleDepthSame(const Data2D &data2D, const Data2D &filter) const;

    private:
        Data4D          m_kernels;         ///< Kernels
        Data1D          m_bias;            ///< Bias
        std::string     m_borderMode;      ///< Method to deal with boardering squares
        int             m_kernelsCount;    ///< Number of kernels
        int             m_nDepth;          ///< Depth of convolution
        int             m_nRows;           ///< Number of rows in convolution
        int             m_nCols;           ///< Number of columns in convolution
    };

    class LayerDense : public Layer
    {
    public:
        LayerDense();

        void LoadWeights(std::ifstream &inputFileStream);

        DataBlock* CalculateOutput(DataBlock* pDataBlock);

        unsigned int GetNInputRows() const;

        unsigned int GetNInputCols() const;

        unsigned int GetOutputUnits() const;

    private:
        Data2D          m_weights;         ///< Weights to apply
        Data1D          m_bias;            ///< Bias to apply
        int             m_inputCount;      ///< Input count
        int             m_nNeurons;        ///< Number of neurons
    };

    /**
     *  @brief  Load the weights for the Keras model
     */
    void LoadWeights(const std::string &inputFileName);

    /**
     *  @brief  Calculate the outcomes based on the data
     *
     *  @param  pDataBlock pointer to data block
     *  @param  outputData1D output data block to fill
     */
    void CalculateOutput(DataBlock *pDataBlock, Data1D &outputData1D);

    unsigned int GetNInputRows() const;

    unsigned int GetNInputCols() const;

    unsigned int GetOutputLength() const;

    /**
     *  @brief  Load a row vector
     *
     *  @param  inputFileStream
     *  @param  nCols in file
     */
    static Data1D Read1DArray(std::ifstream &inputFileStream, const int nCols);

    typedef std::vector<Layer *> Layers;

    int        m_nLayers;           ///< The number of layers in the model
    int        m_nActiveLayers;     ///< The number of layers in the model that are active in the compute output i.e. excluding dropout layers
    Layers     m_layers;            ///< The layers in the model
    bool       m_verbose;           ///< Print stuff maybe
};

} // namespace lar_content

#endif // #ifndef LAR_KERAS_MODEL
