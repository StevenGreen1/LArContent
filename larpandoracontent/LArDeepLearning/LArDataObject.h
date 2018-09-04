/**
 *  @file   larpandoracontent/LArDeepLearning/LArKerasModel.h
 *
 *  @brief  Header file for the lar keras model class.
 *
 *  $Log: $
 */
#ifndef LAR_DATA_OBJECT_H
#define LAR_DATA_OBJECT_H 1

#include "larpandoracontent/LArHelpers/LArMvaHelper.h"

#include "larpandoracontent/LArObjects/LArMvaInterface.h"

#include "Pandora/Algorithm.h"
#include "Pandora/StatusCodes.h"

#include <functional>
#include <map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  Data1D class
 */
class Data1D
{
public:
    /**
     *  @brief  Default constructor
     */
    Data1D();

    /**
     *  @brief  Constructor
     *
     *  @param  elements to hold in class
     */
    Data1D(const pandora::FloatVector &elements);

    /**
     *  @brief  Constructor
     *
     *  @param  nElements number of elements
     *  @param  init value to initialise elements to
     */
    Data1D(const int nElements, const float init);

    /**
     *  @brief  Append element to object vector
     *
     *  @param  element to append
     */
    void Append(const float element);

    /**
     *  @brief  Set element in the object if it exists
     *
     *  @param  idxI index of element to set
     *  @param  element to set element to
     */
    void Set(int idxI, float element);

    /**
     *  @brief  Get element in the object if it exists
     *
     *  @param  idxI index of element to get
     */
    float Get(int idxI) const;

    /**
     *  @brief  Get size of data object
     */
    unsigned int GetSizeI() const;

private:
    pandora::FloatVector     m_elements;     ///< Data
    int                      m_sizeI;      ///< Number of elements
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data1D::Data1D() :
    m_sizeI(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data1D::Data1D(const pandora::FloatVector &elements) :
    m_elements(elements),
    m_sizeI(elements.size())
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data1D::Data1D(const int nElements, const float init) :
    m_elements(pandora::FloatVector(nElements, init)),
    m_sizeI(nElements)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Data1D::Append(float element)
{
    m_elements.push_back(element);
    m_sizeI++;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Data1D::Set(int idxI, float element)
{
    if (idxI < m_sizeI)
    {
        m_elements.at(idxI) = element;
        return;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Data1D::Get(int idxI) const
{
    if (idxI < m_sizeI)
        return m_elements.at(idxI);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data1D::GetSizeI() const
{
    return m_sizeI;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Data2D class
 */
class Data2D
{
public:
    /**
     *  @brief  Default constructor
     */
    Data2D();

    /**
     *  @brief  Constructor
     *
     *  @param  elements to hold in class
     */
    Data2D(const std::vector<Data1D> &elements);

    /**
     *  @brief  Constructor
     *
     *  @param  nElements number of elements
     *  @param  init value to initialise elements to
     */
    Data2D(const int nElements, const Data1D &init);

    /**
     *  @brief  Append element to object vector
     *
     *  @param  element to append
     */
    void Append(const Data1D &element);

    /**
     *  @brief  Set element in the object if it exists
     *
     *  @param  idxI first index of element to set
     *  @param  idxJ second index of element to set
     *  @param  element to set element to
     */
    void Set(int idxI, int idxJ, float element);

    /**
     *  @brief  Get element in the object if it exists
     *
     *  @param  idxI index of element to get
     *  @param  idxJ index of element to get
     */
    float Get(int idxI, int idxJ) const;

    /**
     *  @brief  Get Data1D element in the object if it exists
     *
     *  @param  idxJ index of element to get
     */
    Data1D GetData1D(int idxJ) const;

    /**
     *  @brief  Get size of data object along first dimension
     */
    unsigned int GetSizeI() const;

    /**
     *  @brief  Get size of data object along second dimension
     */
    unsigned int GetSizeJ() const;

private:
    std::vector<Data1D>      m_elements;   ///< Data
    int                      m_sizeI;      ///< Number of elements along first dimension
    int                      m_sizeJ;      ///< Number of elements along second dimension
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data2D::Data2D() :
    m_sizeJ(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data2D::Data2D(const std::vector<Data1D> &elements) :
    m_elements(elements),
    m_sizeJ(elements.size())
{
    if (m_sizeJ > 0)
        m_sizeI = elements.front().GetSizeI();

    for (const Data1D &data1D : elements)
    {
        if (data1D.GetSizeI() != m_sizeI)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data2D::Data2D(const int nElements, const Data1D &init) :
    m_elements(std::vector<Data1D>(nElements, init)),
    m_sizeI(0),
    m_sizeJ(nElements)
{
    if (m_sizeJ > 0)
        m_sizeI = init.GetSizeI();
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Data2D::Append(const Data1D &element)
{
    if (m_sizeJ == 0)
        m_sizeI = element.GetSizeI();

    if (m_sizeI != element.GetSizeI())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    m_elements.push_back(element);
    m_sizeJ++;
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Data2D::Set(int idxI, int idxJ, float element)
{
    if (idxJ < m_sizeJ)
    {
        m_elements.at(idxJ).Set(idxI, element);
        return;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Data2D::Get(int idxI, int idxJ) const
{
    if (idxJ < m_sizeJ)
        return m_elements.at(idxJ).Get(idxI);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data1D Data2D::GetData1D(int idxJ) const
{
    if (idxJ < m_sizeJ)
        return m_elements.at(idxJ);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data2D::GetSizeI() const
{
    return m_sizeI;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data2D::GetSizeJ() const
{
    return m_sizeJ;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Data3D class
 */
class Data3D
{
public:
    /**
     *  @brief  Default constructor
     */
    Data3D();

    /**
     *  @brief  Constructor
     *
     *  @param  elements to hold in class
     */
    Data3D(const std::vector<Data2D> &elements);

    /**
     *  @brief  Constructor
     *
     *  @param  nElements number of elements
     *  @param  init value to initialise elements to
     */
    Data3D(const int nElements, const Data2D &init);

    /**
     *  @brief  Append element to object vector
     *
     *  @param  element to append
     */
    void Append(const Data2D &element);

    /**
     *  @brief  Set element in the object if it exists
     *
     *  @param  idxI first index of element to set
     *  @param  idxJ second index of element to set
     *  @param  idxK third index of element to set
     *  @param  element to set element to
     */
    void Set(int idxI, int idxJ, int idxK, float element);

    /**
     *  @brief  Get element in the object if it exists
     *
     *  @param  idxI index of element to get
     *  @param  idxJ index of element to get
     *  @param  idxK index of element to get
     */
    float Get(int idxI, int idxJ, int idxK) const;

    /**
     *  @brief  Get Data1D element in the object if it exists
     *
     *  @param  idxJ index of element to get
     *  @param  idxK index of element to get
     */
    Data1D GetData1D(int idxJ, int idxK) const;

    /**
     *  @brief  Get Data2D element in the object if it exists
     *
     *  @param  idxK index of element to get
     */
    Data2D GetData2D(int idxK) const;

    /**
     *  @brief  Get size of data object along first dimension
     */
    unsigned int GetSizeI() const;

    /**
     *  @brief  Get size of data object along second dimension
     */
    unsigned int GetSizeJ() const;

    /**
     *  @brief  Get size of data object along third dimension
     */
    unsigned int GetSizeK() const;

private:
    std::vector<Data2D>      m_elements;     ///< Data
    int                      m_sizeI;      ///< Number of elements along first dimension
    int                      m_sizeJ;      ///< Number of elements along second dimension
    int                      m_sizeK;      ///< Number of elements along third dimension
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data3D::Data3D() :
    m_sizeK(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data3D::Data3D(const std::vector<Data2D> &elements) :
    m_elements(elements),
    m_sizeK(elements.size())
{
    if (m_sizeK > 0)
    {
        m_sizeJ = elements.front().GetSizeJ();
        m_sizeI = elements.front().GetSizeI();
    }

    for (const Data2D &data2D : elements)
    {
        if (data2D.GetSizeJ() != m_sizeJ || data2D.GetSizeI() != m_sizeI)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data3D::Data3D(const int nElements, const Data2D &init) :
    m_elements(std::vector<Data2D>(nElements, init)),
    m_sizeI(0),
    m_sizeJ(0),
    m_sizeK(nElements)
{
    if (m_sizeK > 0)
    {
        m_sizeJ = init.GetSizeJ();
        m_sizeI = init.GetSizeI();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Data3D::Append(const Data2D &element)
{
    if (m_sizeK == 0)
    {
        m_sizeJ = element.GetSizeJ();
        m_sizeI = element.GetSizeI();
    }

    if (m_sizeJ != element.GetSizeJ() || m_sizeI != element.GetSizeI())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    m_elements.push_back(element);
    m_sizeK++;
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Data3D::Set(int idxI, int idxJ, int idxK, float element)
{
    if (idxK < m_sizeK)
    {
        m_elements.at(idxK).Set(idxI, idxJ, element);
        return;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Data3D::Get(int idxI, int idxJ, int idxK) const
{
    if (idxK < m_sizeK)
        return m_elements.at(idxK).Get(idxI, idxJ);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data1D Data3D::GetData1D(int idxJ, int idxK) const
{
   if (idxK < m_sizeK)
        return m_elements.at(idxK).GetData1D(idxJ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data2D Data3D::GetData2D(int idxK) const
{
    if (idxK < m_sizeK)
        return m_elements.at(idxK);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data3D::GetSizeI() const
{
    return m_sizeI;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data3D::GetSizeJ() const
{
    return m_sizeJ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data3D::GetSizeK() const
{
    return m_sizeK;
}

//------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  Data4D class
 */
class Data4D
{
public:
    /**
     *  @brief  Default constructor
     */
    Data4D();

    /**
     *  @brief  Constructor
     *
     *  @param  elements to hold in class
     */
    Data4D(const std::vector<Data3D> &elements);

    /**
     *  @brief  Constructor
     *
     *  @param  nElements number of elements
     *  @param  init value to initialise elements to
     */
    Data4D(const int nElements, const Data3D &init);

    /**
     *  @brief  Append element to object vector
     *
     *  @param  element to append
     */
    void Append(const Data3D &element);

    /**
     *  @brief  Set element in the object if it exists
     *
     *  @param  idxI first index of element to set
     *  @param  idxJ second index of element to set
     *  @param  idxK third index of element to set
     *  @param  idxL fourth index of element to set
     *  @param  element to set element to
     */
    void Set(int idxI, int idxJ, int idxK, int idxL, float element);

    /**
     *  @brief  Get element in the object if it exists
     *
     *  @param  idxI index of element to get
     *  @param  idxJ index of element to get
     *  @param  idxK index of element to get
     *  @param  idxL index of element to get
     */
    float Get(int idxI, int idxJ, int idxK, int idxL) const;

    /**
     *  @brief  Get Data2D element in the object if it exists
     *
     *  @param  idxJ index of element to get
     *  @param  idxK index of element to get
     *  @param  idxL index of element to get
     */
    Data1D GetData1D(int idxJ, int idxK, int idxL) const;

    /**
     *  @brief  Get Data2D element in the object if it exists
     *
     *  @param  idxK index of element to get
     *  @param  idxL index of element to get
     */
    Data2D GetData2D(int idxK, int idxL) const;

    /**
     *  @brief  Get Data3D element in the object if it exists
     *
     *  @param  idxL index of element to get
     */
    Data3D GetData3D(int idxL) const;

    /**
     *  @brief  Get size of data object along first dimension
     */
    unsigned int GetSizeI() const;

    /**
     *  @brief  Get size of data object along second dimension
     */
    unsigned int GetSizeJ() const;

    /**
     *  @brief  Get size of data object along third dimension
     */
    unsigned int GetSizeK() const;

    /**
     *  @brief  Get size of data object along fourth dimension
     */
    unsigned int GetSizeL() const;

private:
    std::vector<Data3D>      m_elements;     ///< Data
    int                      m_sizeI;      ///< Number of elements along first dimension
    int                      m_sizeJ;      ///< Number of elements along second dimension
    int                      m_sizeK;      ///< Number of elements along third dimension
    int                      m_sizeL;      ///< Number of elements along third dimension
};

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data4D::Data4D() :
    m_sizeL(0)
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data4D::Data4D(const std::vector<Data3D> &elements) :
    m_elements(elements),
    m_sizeL(elements.size())
{
    if (m_sizeL > 0)
    {
        m_sizeK = elements.front().GetSizeK();
        m_sizeJ = elements.front().GetSizeJ();
        m_sizeI = elements.front().GetSizeI();
    }

    for (const Data3D &data3D : elements)
    {
        if (data3D.GetSizeK() != m_sizeK || data3D.GetSizeJ() != m_sizeJ || data3D.GetSizeI() != m_sizeI)
            throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data4D::Data4D(const int nElements, const Data3D &init) :
    m_elements(std::vector<Data3D>(nElements, init)),
    m_sizeI(0),
    m_sizeJ(0),
    m_sizeK(0),
    m_sizeL(nElements)
{
    if (m_sizeL > 0)
    {
        m_sizeK = init.GetSizeK();
        m_sizeJ = init.GetSizeJ();
        m_sizeI = init.GetSizeI();
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Data4D::Append(const Data3D &element)
{
    if (m_sizeL == 0)
    {
        m_sizeK = element.GetSizeK();
        m_sizeJ = element.GetSizeJ();
        m_sizeI = element.GetSizeI();
    }

    if (m_sizeK != element.GetSizeK() || m_sizeJ != element.GetSizeJ() || m_sizeI != element.GetSizeI())
        throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);

    m_elements.push_back(element);
    m_sizeL++;
    return;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline void Data4D::Set(int idxI, int idxJ, int idxK, int idxL, float element)
{
    if (idxL < m_sizeL)
    {
        m_elements.at(idxL).Set(idxI, idxJ, idxK, element);
        return;
    }

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline float Data4D::Get(int idxI, int idxJ, int idxK, int idxL) const
{
    if (idxL < m_sizeL)
        return m_elements.at(idxL).Get(idxI, idxJ, idxK);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data1D Data4D::GetData1D(int idxJ, int idxK, int idxL) const
{
    if (idxL < m_sizeL)
        return m_elements.at(idxL).GetData1D(idxJ, idxK);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data2D Data4D::GetData2D(int idxK, int idxL) const
{
    if (idxL < m_sizeL)
        return m_elements.at(idxL).GetData2D(idxK);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline Data3D Data4D::GetData3D(int idxL) const
{
    if (idxL < m_sizeL)
        return m_elements.at(idxL);

    throw pandora::StatusCodeException(pandora::STATUS_CODE_INVALID_PARAMETER);
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data4D::GetSizeI() const
{
    return m_sizeI;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data4D::GetSizeJ() const
{
    return m_sizeJ;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data4D::GetSizeK() const
{
    return m_sizeK;
}

//------------------------------------------------------------------------------------------------------------------------------------------

inline unsigned int Data4D::GetSizeL() const
{
    return m_sizeL;
}

//------------------------------------------------------------------------------------------------------------------------------------------


} // namespace lar_content

#endif // #ifndef LAR_DATA_OBJECT_H
