export class Vector2 {
    static get zero(): Vector2;
    static get right(): Vector2;
    static get up(): Vector2;
    static add(vec1: Vector2, vec2: Vector2): Vector2;
    static sub(vec1: Vector2, vec2: Vector2): Vector2;
    static mul(vec1: Vector2, scalar: number): Vector2;
    static div(vec1: Vector2, scalar: number): Vector2;
    /**
     * Returns dot product of two vectors
     *
     * @static
     * @param {Vector2} vec1
     * @param {Vector2} vec2
     * @returns {number}
     * @memberof Vector2
     */
    static dot(vec1: Vector2, vec2: Vector2): number;
    /**
     * Returns cross product of two vectors
     *
     * @static
     * @param {Vector2} vec1
     * @param {Vector2} vec2
     * @returns {number}
     * @memberof Vector2
     */
    static cross(vec1: Vector2, vec2: Vector2): number;
    /**
     * Returns the angle between two vectors in radians
     *
     * @static
     * @param {Vector2} vec1
     * @param {Vector2} vec2
     * @returns
     * @memberof Vector2
     */
    static angleBetween(vec1: Vector2, vec2: Vector2): number;
    /**
     * Returns the number which is the distance between two vectors
     * @param {Vector2} vector1
     * @param {Vector2} vector2
     */
    static distance(vector1: Vector2, vector2: Vector2): number;
    /**
     * Linearly interpolates between two vectors
     *
     * @param {Vector2} vec1
     * @param {Vector2} vec2
     * @param {number} t Interpolation factor - should be in range [0; 1]
     * @return {Vector2}
     */
    static lerp(vec1: Vector2, vec2: Vector2, t: number): Vector2;
    /**
     * @static
     * @param {number[]} array
     * @returns
     * @memberof Vector2
     */
    static fromArray(array: number[]): Vector2;
    /**
     * Creates an instance of Vector3.
     * @param {number} x
     * @param {number} y
     * @memberof Vector2
     */
    constructor(x: number, y: number);
    get x(): number;
    get y(): number;
    get magnitude(): number;
    get magnitudeSquared(): number;
    add(vec: Vector2): Vector2;
    sub(vec: Vector2): Vector2;
    /**
     *
     *
     * @param {number} scalar
     * @memberof Vector2
     */
    mul(scalar: number): Vector2;
    /**
     *
     *
     * @param {number} scalar
     * @memberof Vector2
     */
    div(scalar: number): Vector2;
    /**
     * Returns dot product of two vectors
     *
     * @static
     * @param {Vector2} vec
     * @returns {number}
     * @memberof Vector2
     */
    dot(vec: Vector2): number;
    /**
     * Returns cross product of two vectors
     *
     * @param {Vector2} vec
     * @returns {number}
     * @memberof Vector2
     */
    cross(vec: Vector2): number;
    /**
     * Returns the angle between two vectors in radians
     *
     * @static
     * @param {Vector2} vec
     * @returns
     * @memberof Vector2
     */
    angleBetween(vec: Vector2): number;
    /**
     * Returns the number which is the distance between two vectors
     * @param {Vector2} vector
     */
    distance(vector: Vector2): number;
    /**
     * Returns the vector, the madnitude of which is limited by specified one
     *
     * @param {number} maxLength
     * @returns
     */
    clampMagnitude(maxLength: number): Vector2;
    /**
     * Linearly interpolates between two vectors
     *
     * @param {Vector2} vector
     * @param {number} t Interpolation factor - should be in range [0; 1]
     * @return {Vector2}
     */
    lerp(vector: Vector2, t: number): Vector2;
    neg(): Vector2;
    normalize(): Vector2;
    toPoint2dSignal(): void;
    toString(digits?: number): string;
    /**
     * @return {[number, number]}
     * @memberof Vector2
     */
    toArray(): [number, number];
}
export default Vector2;
