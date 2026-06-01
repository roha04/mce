#include "MathCore.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

// Реалізація MathCore: від сітки HEX20 до локальних K^e та F^e (етапи 1–2 методички).

namespace
{
// ---------------------------------------------------------------------------
// Допоміжні дані HEX20 (не експортуються назовні)
// ---------------------------------------------------------------------------

// Індекси вузлів на «логічній» сітці 0..2Nx по кожній осі (0,1,2 = кут/середина ребра).
// Для localNode: ix = 2*ex + offset[0], аналогічно iy, iz — зв'язок локального № з сіткою.
constexpr int kHex20Offsets[20][3] = {
    {0, 0, 0}, {2, 0, 0}, {2, 2, 0}, {0, 2, 0},
    {0, 0, 2}, {2, 0, 2}, {2, 2, 2}, {0, 2, 2},
    {1, 0, 0}, {2, 1, 0}, {1, 2, 0}, {0, 1, 0},
    {0, 0, 1}, {2, 0, 1}, {2, 2, 1}, {0, 2, 1},
    {1, 0, 2}, {2, 1, 2}, {1, 2, 2}, {0, 1, 2}};

// Тип вузла в локальній системі (ξ,η,ζ) ∈ [-1,1]³ — визначає вид формули ψ_i.
enum class ShapeNodeKind
{
    Corner,   // 8 кутів: трилінійний serendipity
    EdgeXi,   // 4 ребра, паралельні осі ξ (ξ=0 на ребрі)
    EdgeEta,  // 4 ребра, паралельні η
    EdgeZeta  // 4 ребра, паралельні ζ
};

// Метадані для i-го вузла: kind і знаки s0,s1,s2 у (1±s·координата).
struct ShapeNodeMeta
{
    ShapeNodeKind kind;
    int s0;
    int s1;
    int s2;
};

// Таблиця відповідає порядку вузлів у Element::nodes (як у методичці / NT).
constexpr ShapeNodeMeta kShapeNodeMeta[kHex20NodeCount] = {
    {ShapeNodeKind::Corner, -1, -1, -1}, {ShapeNodeKind::Corner, +1, -1, -1},
    {ShapeNodeKind::Corner, +1, +1, -1}, {ShapeNodeKind::Corner, -1, +1, -1},
    {ShapeNodeKind::Corner, -1, -1, +1}, {ShapeNodeKind::Corner, +1, -1, +1},
    {ShapeNodeKind::Corner, +1, +1, +1}, {ShapeNodeKind::Corner, -1, +1, +1},
    {ShapeNodeKind::EdgeXi, 0, -1, -1},  {ShapeNodeKind::EdgeEta, +1, 0, -1},
    {ShapeNodeKind::EdgeXi, 0, +1, -1},  {ShapeNodeKind::EdgeEta, -1, 0, -1},
    {ShapeNodeKind::EdgeZeta, -1, -1, 0}, {ShapeNodeKind::EdgeZeta, +1, -1, 0},
    {ShapeNodeKind::EdgeZeta, +1, +1, 0}, {ShapeNodeKind::EdgeZeta, -1, +1, 0},
    {ShapeNodeKind::EdgeXi, 0, -1, +1},  {ShapeNodeKind::EdgeEta, +1, 0, +1},
    {ShapeNodeKind::EdgeXi, 0, +1, +1},  {ShapeNodeKind::EdgeEta, -1, 0, +1}};

// Обчислення ψ_i(ξ,η,ζ) та ∂ψ_i/∂ξ, ∂ψ_i/∂η, ∂ψ_i/∂ζ для всіх 20 вузлів.
// Кути — трилінійні добутки (1±ξ)(1±η)(1±ζ)(ξ+η+ζ-2)/8; ребра — квадратичні serendipity.
void evaluateHex20Shape(double xi, double eta, double zeta,
                        double N[kHex20NodeCount],
                        double dN[kHex20NodeCount][3])
{
    for (int i = 0; i < kHex20NodeCount; ++i)
    {
        const ShapeNodeMeta& meta = kShapeNodeMeta[i];
        const double s0 = static_cast<double>(meta.s0);
        const double s1 = static_cast<double>(meta.s1);
        const double s2 = static_cast<double>(meta.s2);

        switch (meta.kind)
        {
        case ShapeNodeKind::Corner:
        {
            // ψ = (1/8)(1+s0·ξ)(1+s1·η)(1+s2·ζ)(s0·ξ+s1·η+s2·ζ−2)
            const double a = 1.0 + s0 * xi;
            const double b = 1.0 + s1 * eta;
            const double c = 1.0 + s2 * zeta;
            const double t = s0 * xi + s1 * eta + s2 * zeta - 2.0;

            N[i] = 0.125 * a * b * c * t;
            dN[i][0] = 0.125 * s0 * b * c * (a + t);
            dN[i][1] = 0.125 * s1 * a * c * (b + t);
            dN[i][2] = 0.125 * s2 * a * b * (c + t);
            break;
        }
        case ShapeNodeKind::EdgeXi:
        {
            // ψ = (1/4)(1−ξ²)(1+s1·η)(1+s2·ζ) — вузол на ребрі ξ=0
            const double a = 1.0 - xi * xi;
            const double b = 1.0 + s1 * eta;
            const double c = 1.0 + s2 * zeta;

            N[i] = 0.25 * a * b * c;
            dN[i][0] = -0.5 * xi * b * c;
            dN[i][1] = 0.25 * a * s1 * c;
            dN[i][2] = 0.25 * a * b * s2;
            break;
        }
        case ShapeNodeKind::EdgeEta:
        {
            // ψ = (1/4)(1+s0·ξ)(1−η²)(1+s2·ζ) — вузол на ребрі η=0
            const double a = 1.0 + s0 * xi;
            const double b = 1.0 - eta * eta;
            const double c = 1.0 + s2 * zeta;

            N[i] = 0.25 * a * b * c;
            dN[i][0] = 0.25 * s0 * b * c; // ∂/∂ξ: по (1+s0·ξ), не по (1−η²)
            dN[i][1] = -0.5 * eta * a * c;
            dN[i][2] = 0.25 * a * b * s2;
            break;
        }
        case ShapeNodeKind::EdgeZeta:
        {
            // ψ = (1/4)(1+s0·ξ)(1+s1·η)(1−ζ²) — вузол на ребрі ζ=0
            const double a = 1.0 + s0 * xi;
            const double b = 1.0 + s1 * eta;
            const double c = 1.0 - zeta * zeta;

            N[i] = 0.25 * a * b * c;
            dN[i][0] = 0.25 * b * c * s0;
            dN[i][1] = 0.25 * a * c * s1;
            dN[i][2] = -0.5 * zeta * a * b;
            break;
        }
        }
    }
}

// Матриця деформацій B (6×60): ε = B·U_e, зв'язок градієнтів ψ з компонентами u,v,w.
void buildStrainDisplacementMatrix(const double dNGlobal[kHex20NodeCount][3], double B[6][kElementDofCount])
{
    for (int row = 0; row < 6; ++row)
        for (int col = 0; col < kElementDofCount; ++col)
            B[row][col] = 0.0;

    for (int k = 0; k < kHex20NodeCount; ++k)
    {
        const int base = 3 * k; // DOF вузла k: u=base+0, v=base+1, w=base+2
        const double dNx = dNGlobal[k][0];
        const double dNy = dNGlobal[k][1];
        const double dNz = dNGlobal[k][2];

        // Рядки 0..2: нормальні деформації εxx, εyy, εzz
        B[0][base + 0] = dNx;
        B[1][base + 1] = dNy;
        B[2][base + 2] = dNz;
        // Рядки 3..5: зсувні γxy, γyz, γxz (у Voigt зберігаються як ті самі індекси)
        B[3][base + 0] = dNy;
        B[3][base + 1] = dNx;
        B[4][base + 1] = dNz;
        B[4][base + 2] = dNy;
        B[5][base + 0] = dNz;
        B[5][base + 2] = dNx;
    }
}

// Додає до K^e внесок scale·B^T·D·B з однієї точки Гауса (scale = w·det(J)).
void accumulateBtDB(const double B[6][kElementDofCount],
                    const double D[6][6],
                    double scale,
                    ElementStiffnessMatrix& K)
{
    double DB[6][kElementDofCount]{};

    // Крок 1: DB = D · B (6×60)
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < kElementDofCount; ++j)
        {
            DB[i][j] = 0.0;
            for (int k = 0; k < 6; ++k)
                DB[i][j] += D[i][k] * B[k][j];
        }

    // Крок 2: K += scale · B^T · DB (симетричний внесок однієї точки квадратури)
    for (int i = 0; i < kElementDofCount; ++i)
        for (int j = 0; j < kElementDofCount; ++j)
        {
            double sum = 0.0;
            for (int k = 0; k < 6; ++k)
                sum += B[k][i] * DB[k][j];
            K[i][j] += scale * sum;
        }
}
} // namespace

// =============================================================================
// Клас Mesh — етап 1 методички: геометрія, AKT, NT, ZU, ZP
// =============================================================================

// Скидання всіх масивів сітки перед новою генерацією.
void Mesh::clear()
{
    nodes_.clear();
    elements_.clear();
    akt_.clear();
    nt_.clear();
    zu_.clear();
    zp_.clear();
}

// Лінійний індекс вузла на повній логічній сітці (2·Nx+1) × (2·Ny+1) × (2·Nz+1).
// ix,iy,iz — цілі індекси вузла на «дрібній» сітці перед compact-перенумерацією.
int Mesh::gridNodeIndex(int ix, int iy, int iz, int Nx, int Ny)
{
    const int nx = 2 * Nx + 1;
    const int ny = 2 * Ny + 1;
    return ix + iy * nx + iz * nx * ny;
}

// Генерація структурованої сітки: логічна сітка (2·Nx+1)³, потім compact — лише вузли HEX20.
// Уникаємо «сирітських» вузлів, які не входять у жоден 20-вузловий елемент.
void Mesh::generateRectangularParallelepiped(double Lx, double Ly, double Lz,
                                           int Nx, int Ny, int Nz)
{
    if (Nx < 1 || Ny < 1 || Nz < 1)
        throw std::invalid_argument("Nx, Ny, Nz must be >= 1");

    if (Lx <= 0.0 || Ly <= 0.0 || Lz <= 0.0)
        throw std::invalid_argument("Lx, Ly, Lz must be positive");

    clear();

    lx_ = Lx;
    ly_ = Ly;
    lz_ = Lz;
    nx_ = Nx;
    ny_ = Ny;
    nz_ = Nz;

    const int gridNx = 2 * Nx + 1;
    const int gridNy = 2 * Ny + 1;
    const int gridNz = 2 * Nz + 1;
    const int gridSize = gridNx * gridNy * gridNz;

    // --- Крок 1: для кожного СЕ записуємо 20 логічних індексів і позначаємо використані вузли ---
    std::vector<unsigned char> nodeUsed(static_cast<std::size_t>(gridSize), 0);
    std::vector<Element> tempElements(static_cast<std::size_t>(Nx * Ny * Nz));

    std::size_t elemId = 0;
    for (int ez = 0; ez < Nz; ++ez)
    {
        for (int ey = 0; ey < Ny; ++ey)
        {
            for (int ex = 0; ex < Nx; ++ex)
            {
                Element element{};

                for (int local = 0; local < kHex20NodeCount; ++local)
                {
                    const int ix = 2 * ex + kHex20Offsets[local][0];
                    const int iy = 2 * ey + kHex20Offsets[local][1];
                    const int iz = 2 * ez + kHex20Offsets[local][2];
                    const int logical = gridNodeIndex(ix, iy, iz, Nx, Ny);
                    element.nodes[local] = logical;
                    nodeUsed[static_cast<std::size_t>(logical)] = 1;
                }

                tempElements[elemId++] = element;
            }
        }
    }

    // --- Крок 2: compact — лише вузли, що входять хоча б в один HEX20 (без «сиріт») ---
    std::vector<int> logicalToCompact(static_cast<std::size_t>(gridSize), -1);
    int compactCount = 0;
    for (int iz = 0; iz < gridNz; ++iz)
    {
        for (int iy = 0; iy < gridNy; ++iy)
        {
            for (int ix = 0; ix < gridNx; ++ix)
            {
                const int logical = gridNodeIndex(ix, iy, iz, Nx, Ny);
                if (!nodeUsed[static_cast<std::size_t>(logical)])
                    continue;

                logicalToCompact[static_cast<std::size_t>(logical)] = compactCount++;
            }
        }
    }

    // --- Крок 3: фізичні координати AKT: X = Lx·ix/(2Nx), Y = Ly·iy/(2Ny), Z = Lz·iz/(2Nz) ---
    nodes_.resize(static_cast<std::size_t>(compactCount));
    for (int iz = 0; iz < gridNz; ++iz)
    {
        for (int iy = 0; iy < gridNy; ++iy)
        {
            for (int ix = 0; ix < gridNx; ++ix)
            {
                const int logical = gridNodeIndex(ix, iy, iz, Nx, Ny);
                const int compact = logicalToCompact[static_cast<std::size_t>(logical)];
                if (compact < 0)
                    continue;

                nodes_[static_cast<std::size_t>(compact)] = {
                    Lx * static_cast<double>(ix) / (2.0 * Nx),
                    Ly * static_cast<double>(iy) / (2.0 * Ny),
                    Lz * static_cast<double>(iz) / (2.0 * Nz)};
            }
        }
    }

    // --- Крок 4: перенумерація NT — логічний індекс → компактний глобальний 0..nqp-1 ---
    elements_.resize(tempElements.size());
    for (std::size_t e = 0; e < tempElements.size(); ++e)
    {
        Element element{};
        for (int local = 0; local < kHex20NodeCount; ++local)
        {
            const int logical = tempElements[e].nodes[local];
            element.nodes[local] = logicalToCompact[static_cast<std::size_t>(logical)];
        }
        elements_[e] = element;
    }

    buildAKTandNT();
}

// Формування робочих масивів AKT (7) та NT (21) після побудови nodes_ і elements_.
void Mesh::buildAKTandNT()
{
    const std::size_t nqp = nodes_.size();
    const std::size_t nel = elements_.size();

    // AKT: три блоки по nqp — спочатку всі X, потім Y, потім Z (масив 7 методички).
    akt_.resize(nqp * 3);
    for (std::size_t j = 0; j < nqp; ++j)
    {
        akt_[j] = nodes_[j].x;
        akt_[nqp + j] = nodes_[j].y;
        akt_[2 * nqp + j] = nodes_[j].z;
    }

    // NT: для кожного елемента e — 20 глобальних номерів вузлів (масив 21).
    nt_.resize(nel);
    for (std::size_t e = 0; e < nel; ++e)
    {
        for (int i = 0; i < kHex20NodeCount; ++i)
            nt_[e][static_cast<std::size_t>(i)] = elements_[e].nodes[i];
    }
}

// ZU: вузли нижньої грані Y=0 (закріплення). ZP: верхній шар елементів (ey = ny-1) з тиском P.
void Mesh::buildBoundaryData(const double pressure)
{
    zu_.clear();
    zp_.clear();

    constexpr double kTol = 1.0e-6;

    // ZU (30): закріплена нижня грань — всі вузли з Y ≈ 0 (Ux=Uy=Uz=0 через штраф).
    for (std::size_t nodeId = 0; nodeId < nodes_.size(); ++nodeId)
    {
        if (std::abs(nodes_[nodeId].y) <= kTol)
            zu_.push_back(static_cast<int>(nodeId));
    }

    // ZP (31): тиск P на верхній грані — елементи верхнього шару (ey = Ny−1).
    for (int ez = 0; ez < nz_; ++ez)
    {
        for (int ey = 0; ey < ny_; ++ey)
        {
            for (int ex = 0; ex < nx_; ++ex)
            {
                if (ey != ny_ - 1)
                    continue;

                const std::size_t elemIndex =
                    static_cast<std::size_t>(ex + ey * nx_ + ez * nx_ * ny_);

                ZPEntry entry;
                entry.elementIndex = static_cast<int>(elemIndex);
                entry.faceIndex = kHexFaceEtaPlus;
                entry.pressure = pressure;
                zp_.push_back(entry);
            }
        }
    }
}

// Півширина стрічки ng для глобальної матриці MG (формула 22, заняття 9).
// Обчислюємо максимум (N_max − N_min) по всіх СЕ, множимо на 3 (DOF на вузол).
int Mesh::computeHalfBandwidth() const
{
    if (elements_.empty())
        return 0;

    int maxNodeSpan = 0;

    for (const Element& element : elements_)
    {
        // A_{0,e} — початкові min/max для поточного СЕ.
        int nMin = element.nodes[0];
        int nMax = element.nodes[0];

        // Проходимо всі 20 локальних вузлів і порівнюємо глобальні індекси A_{i,e}.
        for (int i = 1; i < kHex20NodeCount; ++i)
        {
            const int globalNode = element.nodes[i];
            nMin = std::min(nMin, globalNode);
            nMax = std::max(nMax, globalNode);
        }

        maxNodeSpan = std::max(maxNodeSpan, nMax - nMin);
    }

    // 3 DOF на вузол: max |dof_i - dof_j| = 3*(nMax-nMin) + 2 (u..w).
    // Стрічка зберігає пари з j-i < halfBandwidth, тому L >= 3*span + 3.
    return (maxNodeSpan + 1) * 3;
}

// =============================================================================
// Квадратура Гауса в об'ємі елемента (27 точок = 3×3×3)
// =============================================================================
// Вузли: ±√(3/5) та 0; ваги 5/9, 8/9, 5/9 — точна інтеграція поліномів до 5-го степеня.
// Використовується при ∫ B^T D B det(J) dV у assembleElementStiffnessMatrix.
std::vector<GaussPoint3D> generateGaussPoints3D()
{
    static const double coords[3] = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
    static const double weights[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

    std::vector<GaussPoint3D> points;
    points.reserve(kGaussPointCount3D);

    for (int k = 0; k < 3; ++k)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int i = 0; i < 3; ++i)
            {
                GaussPoint3D gp;
                gp.xi = coords[i];
                gp.eta = coords[j];
                gp.zeta = coords[k];
                gp.weight = weights[i] * weights[j] * weights[k];
                points.push_back(gp);
            }
        }
    }

    return points;
}

// =============================================================================
// Базисні функції HEX20 (публічні обгортки)
// =============================================================================

// Повертає лише значення ψ_i(ξ,η,ζ) у масив N[20] (похідні обчислюються всередині, але не повертаються).
void computeHex20ShapeFunctions(double xi, double eta, double zeta, double N[kHex20NodeCount])
{
    double dN[kHex20NodeCount][3];
    evaluateHex20Shape(xi, eta, zeta, N, dN);
}

// Повертає ∂ψ_i/∂ξ, ∂ψ_i/∂η, ∂ψ_i/∂ζ у локальних координатах (для J та B).
void computeHex20ShapeDerivativesLocal(double xi, double eta, double zeta,
                                       double dN[kHex20NodeCount][3])
{
    double N[kHex20NodeCount];
    evaluateHex20Shape(xi, eta, zeta, N, dN);
}

// =============================================================================
// Кеш DFIABG (39) — похідні базису в 27 вузлах Гауса об'єму
// =============================================================================
// Обчислюється один раз при старті програми; не залежить від Lx,Ly,Lz та NT.
GaussShapeDerivativeCache::GaussShapeDerivativeCache()
{
    points_ = generateGaussPoints3D();

    for (int gp = 0; gp < kGaussPointCount3D; ++gp)
    {
        double N[kHex20NodeCount];
        evaluateHex20Shape(points_[static_cast<std::size_t>(gp)].xi,
                           points_[static_cast<std::size_t>(gp)].eta,
                           points_[static_cast<std::size_t>(gp)].zeta,
                           N,
                           dN_[gp]);
    }
}

// Singleton: спільний кеш для всіх елементів (економія часу збірки K^e).
const GaussShapeDerivativeCache& GaussShapeDerivativeCache::instance()
{
    static const GaussShapeDerivativeCache cache;
    return cache;
}

// Вказівник на [∂ψ/∂ξ, ∂ψ/∂η, ∂ψ/∂ζ] для вузла node у точці Гауса gp.
const double* GaussShapeDerivativeCache::dNLocal(int gp, int node) const
{
    return dN_[gp][node];
}

// =============================================================================
// Кеш DPSITE (46) — базис на грані η = +1 (верхня грань у нашій постановці)
// =============================================================================
// 9 точок Гауса 2D; h = ξ, t = ζ; η зафіксовано в +1 (верхня грань паралелепіпеда).
FaceGaussShapeCache::FaceGaussShapeCache()
{
    points2d_ = generateGaussPoints2D();
    constexpr double etaFace = 1.0;

    for (int gp = 0; gp < kGaussPointCount2D; ++gp)
    {
        const GaussPoint3D& point = points2d_[static_cast<std::size_t>(gp)];
        double dN[kHex20NodeCount][3];
        evaluateHex20Shape(point.xi, etaFace, point.eta, psi_[gp], dN);

        for (int node = 0; node < kHex20NodeCount; ++node)
        {
            dPsiDh_[gp][node] = dN[node][0];
            dPsiDt_[gp][node] = dN[node][2];
        }
    }
}

const FaceGaussShapeCache& FaceGaussShapeCache::instance()
{
    static const FaceGaussShapeCache cache;
    return cache;
}

// ψ_i на грані в точці Гауса gp (для згинання тиску у F^e).
double FaceGaussShapeCache::shapeValue(const int gp, const int node) const
{
    return psi_[gp][node];
}

// ∂ψ_i/∂h = ∂ψ_i/∂ξ на грані η=+1.
double FaceGaussShapeCache::shapeDerivH(const int gp, const int node) const
{
    return dPsiDh_[gp][node];
}

// ∂ψ_i/∂t = ∂ψ_i/∂ζ на грані η=+1.
double FaceGaussShapeCache::shapeDerivT(const int gp, const int node) const
{
    return dPsiDt_[gp][node];
}

// =============================================================================
// Якобіан та перехід локальні → глобальні похідні
// =============================================================================

// Обернення 3×3 методом алгебраїчних доповнень; inverse = J^{-1}.
// При |det(J)| < 1e-14 елемент вироджений (перевертений або поганий порядок вузлів NT).
bool Jacobian3x3::invert(double inverse[3][3]) const
{
    const double a = data[0][0];
    const double b = data[0][1];
    const double c = data[0][2];
    const double d = data[1][0];
    const double e = data[1][1];
    const double f = data[1][2];
    const double g = data[2][0];
    const double h = data[2][1];
    const double i = data[2][2];

    const double A = e * i - f * h;
    const double B = -(d * i - f * g);
    const double C = d * h - e * g;
    const double D = -(b * i - c * h);
    const double E = a * i - c * g;
    const double F = -(a * h - b * g);
    const double G = b * f - c * e;
    const double H = -(a * f - c * d);
    const double I = a * e - b * d;

    const double det = a * A + b * B + c * C;
    if (std::abs(det) < 1.0e-14)
        return false;

    const double invDet = 1.0 / det;
    inverse[0][0] = A * invDet;
    inverse[0][1] = D * invDet;
    inverse[0][2] = G * invDet;
    inverse[1][0] = B * invDet;
    inverse[1][1] = E * invDet;
    inverse[1][2] = H * invDet;
    inverse[2][0] = C * invDet;
    inverse[2][1] = F * invDet;
    inverse[2][2] = I * invDet;
    return true;
}

// J_ij = Σ_k x_k,i · ∂ψ_k/∂ξ_j — перехід від локальних до фізичних координат.
Jacobian3x3 computeJacobian(const double dNLocal[kHex20NodeCount][3],
                            const Node elementNodes[kHex20NodeCount])
{
    Jacobian3x3 jacobian{};

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
            jacobian.data[i][j] = 0.0;
    }

    // J[i][j] = ∂x_i/∂ξ_j = Σ_k X_k,i · ∂ψ_k/∂ξ_j  (i,j — рядок/стовпець фізична/локальна)
    for (int k = 0; k < kHex20NodeCount; ++k)
    {
        const double coords[3] = {elementNodes[k].x, elementNodes[k].y, elementNodes[k].z};

        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                jacobian.data[i][j] += coords[i] * dNLocal[k][j];
    }

    // det(J) для множника w·det(J) у інтегралі об'єму
    const double a = jacobian.data[0][0];
    const double b = jacobian.data[0][1];
    const double c = jacobian.data[0][2];
    const double d = jacobian.data[1][0];
    const double e = jacobian.data[1][1];
    const double f = jacobian.data[1][2];
    const double g = jacobian.data[2][0];
    const double h = jacobian.data[2][1];
    const double i = jacobian.data[2][2];

    jacobian.determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
    return jacobian;
}

// ∂ψ/∂x_k = Σ_m (J^{-1})_{km} · ∂ψ/∂ξ_m — перехід до фізичних похідних для матриці B.
void transformShapeDerivativesToGlobal(const double dNLocal[kHex20NodeCount][3],
                                       const double invJ[3][3],
                                       double dNGlobal[kHex20NodeCount][3])
{
    for (int k = 0; k < kHex20NodeCount; ++k)
    {
        const double dNdXi = dNLocal[k][0];
        const double dNdEta = dNLocal[k][1];
        const double dNdZeta = dNLocal[k][2];

        dNGlobal[k][0] = invJ[0][0] * dNdXi + invJ[1][0] * dNdEta + invJ[2][0] * dNdZeta;
        dNGlobal[k][1] = invJ[0][1] * dNdXi + invJ[1][1] * dNdEta + invJ[2][1] * dNdZeta;
        dNGlobal[k][2] = invJ[0][2] * dNdXi + invJ[1][2] * dNdEta + invJ[2][2] * dNdZeta;
    }
}

// D для ізотропного тіла: λ, μ через E та ν (плоский/простірний закон Гука у векторній формі Voigt).
void buildIsotropicElasticityMatrix6x6(const MaterialProperties& material, double D[6][6])
{
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            D[i][j] = 0.0;

    const double E = material.youngModulus;
    const double nu = material.poissonRatio;
    const double factor = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
    const double lambda = factor * nu;  // перша Lamé-константа
    const double mu = E / (2.0 * (1.0 + nu));  // модуль зсуву G

    // Діагональні блоки нормальних напружень; позадіагональ — λ (об'ємний зв'язок)
    D[0][0] = D[1][1] = D[2][2] = lambda + 2.0 * mu;
    D[0][1] = D[0][2] = D[1][0] = D[1][2] = D[2][0] = D[2][1] = lambda;
    D[3][3] = D[4][4] = D[5][5] = mu;
}

// =============================================================================
// Збірка локальної матриці жорсткості K^e (60×60)
// =============================================================================
// K^e = Σ_gp w_gp · det(J) · B^T D B; перевірка det(J)>0; симетризація після інтегрування.
ElementStiffnessMatrix assembleElementStiffnessMatrix(
    const Node elementNodes[kHex20NodeCount],
    const MaterialProperties& material,
    const GaussShapeDerivativeCache& gaussCache)
{
    ElementStiffnessMatrix K{};
    for (auto& row : K)
        row.fill(0.0);

    double D[6][6];
    buildIsotropicElasticityMatrix6x6(material, D);

    const auto& points = gaussCache.gaussPoints();

    for (int gp = 0; gp < kGaussPointCount3D; ++gp)
    {
        // 1) Взяти ∂ψ/∂(ξ,η,ζ) з DFIABG (не перераховувати щоразу)
        double dNLocal[kHex20NodeCount][3];
        for (int node = 0; node < kHex20NodeCount; ++node)
        {
            const double* localDerivs = gaussCache.dNLocal(gp, node);
            dNLocal[node][0] = localDerivs[0];
            dNLocal[node][1] = localDerivs[1];
            dNLocal[node][2] = localDerivs[2];
        }

        // 2) J та det(J); det(J)≤0 — елемент «перевернутий» або помилка NT
        Jacobian3x3 jacobian = computeJacobian(dNLocal, elementNodes);
        if (jacobian.determinant <= 0.0)
        {
            std::cerr << "ERROR: det(J) = " << jacobian.determinant
                      << " at Gauss point " << gp
                      << " (xi=" << points[static_cast<std::size_t>(gp)].xi
                      << ", eta=" << points[static_cast<std::size_t>(gp)].eta
                      << ", zeta=" << points[static_cast<std::size_t>(gp)].zeta
                      << "). Element may be inverted or node order is wrong.\n";
            throw std::runtime_error("Non-positive Jacobian determinant in element stiffness assembly");
        }

        double invJ[3][3];
        if (!jacobian.invert(invJ))
            throw std::runtime_error("Singular Jacobian in element stiffness assembly");

        // 3) ∂ψ/∂(x,y,z), матриця B, внесок w·det(J)·B^T D B
        double dNGlobal[kHex20NodeCount][3];
        transformShapeDerivativesToGlobal(dNLocal, invJ, dNGlobal);

        double B[6][kElementDofCount];
        buildStrainDisplacementMatrix(dNGlobal, B);

        const double scale = points[static_cast<std::size_t>(gp)].weight * jacobian.determinant;
        accumulateBtDB(B, D, scale, K);
    }

    // Симетризація через усереднення з транспонованою частиною.
    for (int i = 0; i < kElementDofCount; ++i)
        for (int j = i + 1; j < kElementDofCount; ++j)
        {
            const double avg = 0.5 * (K[i][j] + K[j][i]);
            K[i][j] = avg;
            K[j][i] = avg;
        }

    return K;
}

// =============================================================================
// Квадратура Гауса на грані (9 точок = 3×3) — для поверхневого навантаження F^e
// =============================================================================
// Ті самі вузли ±√(3/5) та 0 на квадраті [-1,1]².
// У структурі GaussPoint3D: xi ↔ координата h (вздовж ξ), eta ↔ t (вздовж ζ), zeta не використовується.
// Вага = w_h · w_t — добуток одномірних ваг (інтеграл ∬ ... dh dt).
// Застосовується в assembleElementLoadVectorTopFace разом із кешем DPSITE.
std::vector<GaussPoint3D> generateGaussPoints2D()
{
    static const double coords[3] = {-std::sqrt(3.0 / 5.0), 0.0, std::sqrt(3.0 / 5.0)};
    static const double weights[3] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};

    std::vector<GaussPoint3D> points;
    points.reserve(kGaussPointCount2D);

    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 3; ++i)
        {
            GaussPoint3D gp;
            gp.xi = coords[i];
            gp.eta = coords[j];
            gp.zeta = 0.0;
            gp.weight = weights[i] * weights[j];
            points.push_back(gp);
        }
    }

    return points;
}

// Допоміжно: чи лежить елемент на верхній грані тіла (max Y вузлів ≈ Ly).
bool isElementOnTopFace(const Node elementNodes[kHex20NodeCount], const double ly, const double tolerance)
{
    double maxY = elementNodes[0].y;
    for (int i = 1; i < kHex20NodeCount; ++i)
        maxY = std::max(maxY, elementNodes[i].y);

    return std::abs(maxY - ly) <= tolerance;
}

// Допоміжно: чи лежить елемент на нижній грані (min Y ≈ 0) — для перевірок/фільтрів.
bool isElementOnBottomFace(const Node elementNodes[kHex20NodeCount], const double tolerance)
{
    double minY = elementNodes[0].y;
    for (int i = 1; i < kHex20NodeCount; ++i)
        minY = std::min(minY, elementNodes[i].y);

    return std::abs(minY) <= tolerance;
}

// =============================================================================
// Збірка локального вектора навантаження F^e від тиску на грані
// =============================================================================
// F^e_y = ∫_S ψ_i · t_y dS; t_y = tractionY (наприклад -P); площа через |∂r/∂h × ∂r/∂t|.
ElementLoadVector assembleElementLoadVectorTopFace(const Node elementNodes[kHex20NodeCount],
                                                   const double tractionY,
                                                   const int faceIndex)
{
    ElementLoadVector Fe{};
    Fe.fill(0.0);

    // Підтримується лише грань η=+1 (верх); інші faceIndex — нульовий F^e.
    if (faceIndex != kHexFaceEtaPlus)
        return Fe;

    const FaceGaussShapeCache& dpsite = FaceGaussShapeCache::instance();
    const std::vector<GaussPoint3D>& points = dpsite.gaussPoints2D();

    for (int gp = 0; gp < kGaussPointCount2D; ++gp)
    {
        // Ковзні вектори ∂r/∂h та ∂r/∂t у фізичних координатах (сума по 20 вузлах)
        double tangentH[3]{};
        double tangentT[3]{};
        for (int k = 0; k < kHex20NodeCount; ++k)
        {
            const double coords[3] = {elementNodes[k].x, elementNodes[k].y, elementNodes[k].z};
            const double dH = dpsite.shapeDerivH(gp, k);
            const double dT = dpsite.shapeDerivT(gp, k);
            for (int dim = 0; dim < 3; ++dim)
            {
                tangentH[dim] += dH * coords[dim];
                tangentT[dim] += dT * coords[dim];
            }
        }

        // |∂r/∂h × ∂r/∂t| — елемент площі dS у формулі ∫_S ψ·t dS
        const double normal[3] = {
            tangentH[1] * tangentT[2] - tangentH[2] * tangentT[1],
            tangentH[2] * tangentT[0] - tangentH[0] * tangentT[2],
            tangentH[0] * tangentT[1] - tangentH[1] * tangentT[0]};

        const double areaScale = std::sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
        const double integrationWeight = points[static_cast<std::size_t>(gp)].weight * areaScale;

        // Навантаження лише по Y (компонента v): F_{3i+1} += ψ_i · t_y · w · |...|
        // tractionY зазвичай -P (тиск зверху стискає, тобто негативний напрям Y).
        for (int i = 0; i < kHex20NodeCount; ++i)
            Fe[3 * i + 1] += dpsite.shapeValue(gp, i) * tractionY * integrationWeight;
    }

    return Fe;
}

// =============================================================================
// Перевірка заняття 13: det(J) на еталонному кубі
// =============================================================================
double verifyUnitCubeJacobianDeterminant()
{
    // det(J) = V_фіз / V_станд; V_станд = 8 для куба в [-1,1]³.
    // Для фізичного куба 4×4×4 (вузли в ±2) ізотропна довжина ребра 4 => det(J) ≈ 8 у центрі.
    constexpr double kPhysicalHalfSpan = 2.0;

    Node unitNodes[kHex20NodeCount];
    for (int i = 0; i < kHex20NodeCount; ++i)
    {
        double xi = 0.0;
        double eta = 0.0;
        double zeta = 0.0;
        getHex20LocalNodeCoordinates(i, xi, eta, zeta);
        unitNodes[i].x = kPhysicalHalfSpan * xi;
        unitNodes[i].y = kPhysicalHalfSpan * eta;
        unitNodes[i].z = kPhysicalHalfSpan * zeta;
    }

    double dNLocal[kHex20NodeCount][3];
    computeHex20ShapeDerivativesLocal(0.0, 0.0, 0.0, dNLocal);
    const Jacobian3x3 jacobian = computeJacobian(dNLocal, unitNodes);
    return jacobian.determinant;
}

// Локальні координати вузла localNode: offset 0,1,2 на сітці → ξ,η,ζ ∈ {-1,0,1}.
// Використовується при постпроцесорі напружень (обчислення B у вузлі, а не в точках Гауса).
void getHex20LocalNodeCoordinates(const int localNode, double& xi, double& eta, double& zeta)
{
    xi = static_cast<double>(kHex20Offsets[localNode][0]) - 1.0;
    eta = static_cast<double>(kHex20Offsets[localNode][1]) - 1.0;
    zeta = static_cast<double>(kHex20Offsets[localNode][2]) - 1.0;
}
