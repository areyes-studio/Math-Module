export class Quaternion {
    static get zero(): Quaternion;
    static get one(): Quaternion;
    static get i(): Quaternion;
    static get j(): Quaternion;
    static get k(): Quaternion;
    static get epsilon(): number;
    /**
     *
     *
     * @static
     * @param {number[]} axis
     * @param {number} angle
     * @return {Quaternion}
     * @memberof Quaternion
     */
    static fromAxisAngle(axis: number[], angle: number): Quaternion;
    /**
     *
     * @static
     * @param {Array<number>} u
     * @param {Array<number>} v
     * @return {*}
     * @memberof Quaternion
     */
    static fromBetweenVectors(u: Array<number>, v: Array<number>): any;
    /**
     * @static
     * @return {*}
     * @memberof Quaternion
     */
    static random(): any;
    /**
     * @static
     * @param {number} phi
     * @param {number} theta
     * @param {number} psi
     * @param {string=} order
     * @memberof Quaternion
     */
    static fromEuler(phi: number, theta: number, psi: number, order?: string | undefined): Quaternion;
    /**
     * @static
     * @param {Array<*>} matrix
     * @memberof Quaternion
     */
    static fromMatrix(matrix: Array<any>): Quaternion;
    /**
     * Adds two quaternions Q1 and Q2
     *
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*}
     * @memberof Quaternion
     */
    static add(quat1: Quaternion, quat2: Quaternion): any;
    /**
     * Subtracts two quaternions Q1 and Q2
     *
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*}
     * @memberof Quaternion
     */
    static sub(quat1: Quaternion, quat2: Quaternion): any;
    /**
     * Calculates the additive inverse, or simply it negates the quaternion
     *
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static neg(quat: Quaternion): any;
    /**
     * Calculates the length/modulus/magnitude or the norm of a quaternion
     *
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static norm(quat: Quaternion): any;
    /**
  * Calculates the squared length/modulus/magnitude or the norm of a quaternion
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static normSq(quat: Quaternion): any;
    /**
     * Normalizes the quaternion to have |Q| = 1 as long as the norm is not zero
     * Alternative names are the signum, unit or versor
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static normalize(quat: Quaternion): any;
    /**
     * Calculates the Hamilton product of two quaternions
     * Leaving out the imaginary part results in just scaling the quat
     * @static
     * @param {*} quat1
     * @param {*} quat2
     * @return {*}
     * @memberof Quaternion
     */
    static mul(quat1: Quaternion, quat2: Quaternion): any;
    /**
     * Scales a quaternion by a scalar, faster than using multiplication
     *
     * @static
     * @param {Quaternion} quat
     * @param {number} s
     * @return {*}
     * @memberof Quaternion
     */
    static scale(quat: Quaternion, s: number): any;
    /**
     * Calculates the dot product of two quaternions
     *
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*}
     * @memberof Quaternion
     */
    static dot(quat1: Quaternion, quat2: Quaternion): any;
    /**
     * Calculates the inverse of a quat for non-normalized quats such that
     * Q^-1 * Q = 1 and Q * Q^-1 = 1
     *
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static inverse(quat: Quaternion): any;
    /**
     * Multiplies a quaternion with the inverse of a second quaternion
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*}
     * @memberof Quaternion
     */
    static div(quat1: Quaternion, quat2: Quaternion): any;
    /**
     * Calculates the conjugate of a quaternion
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static conjugate(quat: Quaternion): any;
    /**
     *
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static exp(quat: Quaternion): any;
    /**
     * Calculates the natural logarithm of the quaternion
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static log(quat: Quaternion): any;
    /**
     * Calculates the power of a quaternion raised to a real number or another quaternion
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*}
     * @memberof Quaternion
     */
    static pow(quat1: Quaternion, quat2: Quaternion): any;
    /**
     * Checks if two quats are the same
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*}
     * @memberof Quaternion
     */
    static equals(quat1: Quaternion, quat2: Quaternion): any;
    /**
     * Checks if all parts of a quaternion are finite
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static isFinite(quat: Quaternion): any;
    /**
     * Checks if any of the parts of the quaternion is not a number
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static isNaN(quat: Quaternion): any;
    /**
     * Gets the Quaternion as a well formatted string
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static toString(quat: Quaternion): any;
    /**
     * Returns the real part of the quaternion
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static real(quat: Quaternion): any;
    /**
     * Returns the imaginary part of the quaternion as a 3D vector / array
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static imag(quat: Quaternion): any;
    /**
     * Gets the actual quaternion as a 4D vector / array
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static toVector(quat: Quaternion): any;
    /**
     * Calculates the 3x3 rotation matrix for the current quat
     * @static
     * @param {Quaternion} quat
     * @param {boolean} twoD
     * @return {*}
     * @memberof Quaternion
     */
    static toMatrix(quat: Quaternion, twoD: boolean): any;
    /**
     * Calculates the homogeneous 4x4 rotation matrix for the current quat
     * @static
     * @param {Quaternion} quat
     * @param {boolean} twoD
     * @return {*}
     * @memberof Quaternion
     */
    static toMatrix4(quat: Quaternion, twoD: boolean): any;
    /**
     * Calculates the Euler angles represented by the current quat
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static toEuler(quat: Quaternion): any;
    /**
     * Clones the actual object
     * @static
     * @param {Quaternion} quat
     * @return {*}
     * @memberof Quaternion
     */
    static clone(quat: Quaternion): any;
    /**
     * Rotates a vector according to the current quaternion, assumes |q|=1
     * @link https://www.xarg.org/proof/vector-rotation-using-quaternions/
     * @static
     * @param {Quaternion} quat
     * @param {*} v
     * @return {*}
     * @memberof Quaternion
     */
    static rotateVector(quat: Quaternion, v: any): any;
    /**
     * Gets a function to spherically interpolate between two quaternions
     * @static
     * @param {Quaternion} quat1
     * @param {Quaternion} quat2
     * @return {*}
     * @memberof Quaternion
     */
    static slerp(quat1: Quaternion, quat2: Quaternion): any;
    /**
     * Creates an instance of Vector3.
     * @param {number} w
     * @param {number} x
     * @param {number} y
     * @param {number} z
     * @memberof Vector3
     */
    constructor(w: number, x: number, y: number, z: number);
    get w(): number;
    get x(): number;
    get y(): number;
    get z(): number;
    add(quat: Quaternion): any;
    sub(quat: Quaternion): any;
    neg(): any;
    norm(): any;
    normSq(): any;
    normalize(): any;
    mul(quat: Quaternion): any;
    scale(s: number): any;
    dot(quat: Quaternion): any;
    inverse(): any;
    div(quat: Quaternion): any;
    conjugate(): any;
    exp(): any;
    log(): any;
    pow(quat: Quaternion): any;
    equals(quat: Quaternion): any;
    isFinite(): any;
    isNaN(): any;
    toString(): any;
    real(): any;
    imag(): any;
    toVector(): any;
    toMatrix(twoD: boolean): any;
    toMatrix4(twoD: boolean): any;
    toEuler(): any;
    clone(): any;
    rotateVector(v: []): any;
    slerp(quat: Quaternion): any;
}
export default Quaternion;
