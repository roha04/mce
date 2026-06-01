#pragma once

// =============================================================================
// MathCore — етапи 1–2 методички: сітка, HEX20, інтегрування Гауса, K^e, F^e
// =============================================================================
// Фізична модель: пружний ізотропний паралелепіпед, статика, МСЕ Гальоркіна.
// Тут формуються локальні матриці жорсткості (60×60) та вектори навантаження елемента,
// які далі збираються в глобальну систему K·U = F (див. GlobalSystem).
// =============================================================================

#include <array>
#include <cstddef>
#include <vector>

constexpr int kHex20NodeCount = 20;   // вузлів на один СЕ (serendipity)
constexpr int kElementDofCount = 60;  // 20 вузлів × 3 компоненти переміщення (u,v,w)
constexpr int kGaussPointCount3D = 27; // 3×3×3 точок Гауса для об'єму елемента
constexpr int kGaussPointCount2D = 9;  // 3×3 точок для навантаженої грані

// Номер грані η = +1 у локальних координатах (відповідає верхній грані Y = Ly).
constexpr int kHexFaceEtaPlus = 5;

/// Глобальні координати вузла сітки (м).
struct Node
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/// Один 20-вузловий HEX20: nodes[i] — глобальний номер i-го локального вузла (масив NT).
struct Element
{
    int nodes[kHex20NodeCount]{};
};

/// Параметри ізотропного матеріалу (закон Гука).
struct MaterialProperties
{
    double youngModulus = 2.0e11;  // E, Па
    double poissonRatio = 0.3;   // ν
};

/// Одна точка квадратури Гауса в (ξ,η,ζ) або (ξ,ζ) на грані.
struct GaussPoint3D
{
    double xi = 0.0;
    double eta = 0.0;
    double zeta = 0.0;
    double weight = 0.0;
};

using ElementStiffnessMatrix = std::array<std::array<double, kElementDofCount>, kElementDofCount>;
using ElementLoadVector = std::array<double, kElementDofCount>;

/// Сітка СЕ та робочі масиви методички AKT, NT, ZU, ZP.
class Mesh
{
public:
    /// Запис навантаження ZP (31): який елемент, яка грань, який тиск P [Па].
    struct ZPEntry
    {
        int elementIndex = 0;
        int faceIndex = kHexFaceEtaPlus;
        double pressure = 0.0;
    };

    /// Будує структуровану сітку HEX20 у паралелепіпеді [0,Lx]×[0,Ly]×[0,Lz].
    /// Nx,Ny,Nz — кількість елементів по осях; вузли лише ті, що входять у serendipity-елемент.
    void generateRectangularParallelepiped(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz);

    /// Після generate: заповнює ZU (вузли Y=0) та ZP (верхні грані з тиском P).
    void buildBoundaryData(double pressure);

    /// Півширина стрічки ng за формулою (22) методички: (N_max - N_min + 1) × 3.
    int computeHalfBandwidth() const;

    const std::vector<Node>& getNodes() const { return nodes_; }
    const std::vector<Element>& getElements() const { return elements_; }

    /// AKT (7): координати [X_j, Y_j, Z_j] у плоскому масиві довжини 3·nqp.
    const std::vector<double>& getAKT() const { return akt_; }
    /// NT (21): для елемента e — глобальні номери 20 локальних вузлів.
    const std::vector<std::array<int, kHex20NodeCount>>& getNT() const { return nt_; }
    /// ZU (30): глобальні номери вузлів на закріпленій нижній грані.
    const std::vector<int>& getZU() const { return zu_; }
    /// ZP (31): список навантажених граней.
    const std::vector<ZPEntry>& getZP() const { return zp_; }

    std::size_t nodeCount() const { return nodes_.size(); }
    std::size_t elementCount() const { return elements_.size(); }

    double getLx() const { return lx_; }
    double getLy() const { return ly_; }
    double getLz() const { return lz_; }
    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    int getNz() const { return nz_; }

    void clear();

private:
    static int gridNodeIndex(int ix, int iy, int iz, int Nx, int Ny);
    void buildAKTandNT();

    std::vector<Node> nodes_;
    std::vector<Element> elements_;
    std::vector<double> akt_;
    std::vector<std::array<int, kHex20NodeCount>> nt_;
    std::vector<int> zu_;
    std::vector<ZPEntry> zp_;
    double lx_ = 0.0;
    double ly_ = 0.0;
    double lz_ = 0.0;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
};

/// 27 точок Гауса на кубі [-1,1]³ (ваги 5/9, 8/9, 5/9).
std::vector<GaussPoint3D> generateGaussPoints3D();
std::vector<GaussPoint3D> generateGaussPoints2D();

/// Базисні функції ψ_i та їх похідні в локальних координатах (ξ,η,ζ).
void computeHex20ShapeFunctions(double xi, double eta, double zeta, double N[kHex20NodeCount]);
void computeHex20ShapeDerivativesLocal(double xi, double eta, double zeta,
                                       double dN[kHex20NodeCount][3]);

/// DFIABG (39): таблиця ∂ψ_i/∂(ξ,η,ζ) у 27 вузлах Гауса (не залежить від геометрії).
class GaussShapeDerivativeCache
{
public:
    static const GaussShapeDerivativeCache& instance();

    const std::vector<GaussPoint3D>& gaussPoints() const { return points_; }
    const double* dNLocal(int gp, int node) const;

private:
    GaussShapeDerivativeCache();

    std::vector<GaussPoint3D> points_;
    double dN_[kGaussPointCount3D][kHex20NodeCount][3]{};
};

using DFIABGCache = GaussShapeDerivativeCache;

/// DPSITE (46): ψ та похідні на грані η=+1 у 9 точках Гауса (кеш для F^e).
class FaceGaussShapeCache
{
public:
    static const FaceGaussShapeCache& instance();

    const std::vector<GaussPoint3D>& gaussPoints2D() const { return points2d_; }
    double shapeValue(int gp, int node) const;
    double shapeDerivH(int gp, int node) const;
    double shapeDerivT(int gp, int node) const;

private:
    FaceGaussShapeCache();

    std::vector<GaussPoint3D> points2d_;
    double psi_[kGaussPointCount2D][kHex20NodeCount]{};
    double dPsiDh_[kGaussPointCount2D][kHex20NodeCount]{};
    double dPsiDt_[kGaussPointCount2D][kHex20NodeCount]{};
};

using DPSITECache = FaceGaussShapeCache;

/// Якобіан переходу (ξ,η,ζ) → (x,y,z) та det(J) для інтеграла з вагою.
struct Jacobian3x3
{
    double data[3][3]{};
    double determinant = 0.0;

    bool invert(double inverse[3][3]) const;
};

/// J_ij = ∂x_i/∂ξ_j за координатами вузлів і похідними ψ у локальних координатах.
Jacobian3x3 computeJacobian(const double dNLocal[kHex20NodeCount][3],
                            const Node elementNodes[kHex20NodeCount]);

/// ∂ψ/∂x, ∂ψ/∂y, ∂ψ/∂z через J^{-1}.
void transformShapeDerivativesToGlobal(const double dNLocal[kHex20NodeCount][3],
                                       const double invJ[3][3],
                                       double dNGlobal[kHex20NodeCount][3]);

/// Матриця пружності D (6×6) для вектора деформацій [εxx,εyy,εzz,γxy,γyz,γxz].
void buildIsotropicElasticityMatrix6x6(const MaterialProperties& material, double D[6][6]);

/// Збирає K^e = ∫ B^T D B det(J) dV чисельно (27 точок Гауса), симетризує результат.
ElementStiffnessMatrix assembleElementStiffnessMatrix(
    const Node elementNodes[kHex20NodeCount],
    const MaterialProperties& material,
    const GaussShapeDerivativeCache& gaussCache = GaussShapeDerivativeCache::instance());

/// Збирає F^e на грані η=+1: тяга tractionY [Па] у напрямку Y (зазвичай -P для стиску зверху).
ElementLoadVector assembleElementLoadVectorTopFace(const Node elementNodes[kHex20NodeCount],
                                                   double tractionY,
                                                   int faceIndex = kHexFaceEtaPlus);

bool isElementOnTopFace(const Node elementNodes[kHex20NodeCount], double ly, double tolerance = 1.0e-9);
bool isElementOnBottomFace(const Node elementNodes[kHex20NodeCount], double tolerance = 1.0e-9);

/// Локальні координати (ξ,η,ζ) вузла localNode у [-1,1]³.
void getHex20LocalNodeCoordinates(int localNode, double& xi, double& eta, double& zeta);

/// Перевірка зан. 13: det(J)≈8 у центрі еталонного куба 4×4×4 м.
double verifyUnitCubeJacobianDeterminant();
