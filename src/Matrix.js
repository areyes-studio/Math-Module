import Quaternion from "./Quaternion";
import Vector3 from "./Vector3.js";

export class Matrix {
    constructor(elements = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]) {
        this.elements = elements;
    }

    /**
     * @param {Matrix} matrix
     * @memberof Matrix
     */
    multiply(matrix) {
        // Multiply this matrix by the given matrix
        const elements = [];
        for (let i = 0; i < 16; i++) {
            elements[i] = 0;
            for (let j = 0; j < 4; j++) {
                elements[i] += this.elements[j * 4 + i % 4] * matrix.elements[j + Math.floor(i / 4) * 4];
            }
        }

        return new Matrix(elements)
    }

    /**
     * @param {Vector3} point
     * @return {Vector3} 
     * @memberof Matrix
     */
    transformPoint(point) {
        const elements = this.elements;

        const transformedPoint = Vector3.fromArray([
            point.x * elements[0] + point.y * elements[4] + point.z * elements[8] + elements[12],
            point.x * elements[1] + point.y * elements[5] + point.z * elements[9] + elements[13],
            point.x * elements[2] + point.y * elements[6] + point.z * elements[10] + elements[14],
        ]);

        return transformedPoint;
    }

    /**
     * @param {Quaternion} quaternion
     * @return {Quaternion} 
     * @memberof Matrix
     */
    transformQuaternion(quaternion) {
        const x = quaternion.x;
        const y = quaternion.y;
        const z = quaternion.z;
        const w = quaternion.w;

        const elements = this.elements;

        const array = [
            x * elements[0] + y * elements[4] + z * elements[8] + w * elements[12],
            x * elements[1] + y * elements[5] + z * elements[9] + w * elements[13],
            x * elements[2] + y * elements[6] + z * elements[10] + w * elements[14],
            x * elements[3] + y * elements[7] + z * elements[11] + w * elements[15]
        ]

        const transformedQuaternion = new Quaternion(
            array[3],
            array[0],
            array[1],
            array[2],
        );

        return transformedQuaternion;
    }

    /**
     * @param {Vector3} vector
     * @return {Vector3} 
     * @memberof Matrix
     */
    transformVector(vector) {
        const x = vector.x;
        const y = vector.y;
        const z = vector.z;

        const elements = this.elements;

        const transformedVector = Vector3.fromArray([
            x * elements[0] + y * elements[4] + z * elements[8],
            x * elements[1] + y * elements[5] + z * elements[9],
            x * elements[2] + y * elements[6] + z * elements[10]
        ]);

        return transformedVector;
    }

    inverse() {
        const te = this.elements
        const n11 = te[0], n21 = te[1], n31 = te[2], n41 = te[3];
        const n12 = te[4], n22 = te[5], n32 = te[6], n42 = te[7];
        const n13 = te[8], n23 = te[9], n33 = te[10], n43 = te[11];
        const n14 = te[12], n24 = te[13], n34 = te[14], n44 = te[15];

        const t11 = n23 * n34 * n42 - n24 * n33 * n42 + n24 * n32 * n43 - n22 * n34 * n43 - n23 * n32 * n44 + n22 * n33 * n44;
        const t12 = n14 * n33 * n42 - n13 * n34 * n42 - n14 * n32 * n43 + n12 * n34 * n43 + n13 * n32 * n44 - n12 * n33 * n44;
        const t13 = n13 * n24 * n42 - n14 * n23 * n42 + n14 * n22 * n43 - n12 * n24 * n43 - n13 * n22 * n44 + n12 * n23 * n44;
        const t14 = n14 * n23 * n32 - n13 * n24 * n32 - n14 * n22 * n33 + n12 * n24 * n33 + n13 * n22 * n34 - n12 * n23 * n34;

        const det = n11 * t11 + n21 * t12 + n31 * t13 + n41 * t14;

        if (det === 0) throw new Error('Det === 0')

        const detInv = 1 / det;

        te[0] = t11 * detInv;
        te[1] = (n24 * n33 * n41 - n23 * n34 * n41 - n24 * n31 * n43 + n21 * n34 * n43 + n23 * n31 * n44 - n21 * n33 * n44) * detInv;
        te[2] = (n22 * n34 * n41 - n24 * n32 * n41 + n24 * n31 * n42 - n21 * n34 * n42 - n22 * n31 * n44 + n21 * n32 * n44) * detInv;
        te[3] = (n23 * n32 * n41 - n22 * n33 * n41 - n23 * n31 * n42 + n21 * n33 * n42 + n22 * n31 * n43 - n21 * n32 * n43) * detInv;

        te[4] = t12 * detInv;
        te[5] = (n13 * n34 * n41 - n14 * n33 * n41 + n14 * n31 * n43 - n11 * n34 * n43 - n13 * n31 * n44 + n11 * n33 * n44) * detInv;
        te[6] = (n14 * n32 * n41 - n12 * n34 * n41 - n14 * n31 * n42 + n11 * n34 * n42 + n12 * n31 * n44 - n11 * n32 * n44) * detInv;
        te[7] = (n12 * n33 * n41 - n13 * n32 * n41 + n13 * n31 * n42 - n11 * n33 * n42 - n12 * n31 * n43 + n11 * n32 * n43) * detInv;

        te[8] = t13 * detInv;
        te[9] = (n14 * n23 * n41 - n13 * n24 * n41 - n14 * n21 * n43 + n11 * n24 * n43 + n13 * n21 * n44 - n11 * n23 * n44) * detInv;
        te[10] = (n12 * n24 * n41 - n14 * n22 * n41 + n14 * n21 * n42 - n11 * n24 * n42 - n12 * n21 * n44 + n11 * n22 * n44) * detInv;
        te[11] = (n13 * n22 * n41 - n12 * n23 * n41 - n13 * n21 * n42 + n11 * n23 * n42 + n12 * n21 * n43 - n11 * n22 * n43) * detInv;

        te[12] = t14 * detInv;
        te[13] = (n13 * n24 * n31 - n14 * n23 * n31 + n14 * n21 * n33 - n11 * n24 * n33 - n13 * n21 * n34 + n11 * n23 * n34) * detInv;
        te[14] = (n14 * n22 * n31 - n12 * n24 * n31 - n14 * n21 * n32 + n11 * n24 * n32 + n12 * n21 * n34 - n11 * n22 * n34) * detInv;
        te[15] = (n12 * n23 * n31 - n13 * n22 * n31 + n13 * n21 * n32 - n11 * n23 * n32 - n12 * n21 * n33 + n11 * n22 * n33) * detInv;

        return new Matrix(te);

    }

    /**
     * @param {Vector3} target
     * @param {Vector3} [up=new Vector3(0, 1, 0)]
     * @return {Matrix} 
     * @memberof Matrix
     */
    lookAt(target, up = new Vector3(0, 1, 0)) {
        let z = new Vector3(this.elements[12], this.elements[13], this.elements[14]).sub(target).normalize();
        let x = up.cross(z).normalize();
        let y = z.cross(x);

        const te = this.elements;

        const scale = new Vector3(
            Math.sqrt(te[0] * te[0] + te[1] * te[1] + te[2] * te[2]),
            Math.sqrt(te[4] * te[4] + te[5] * te[5] + te[6] * te[6]),
            Math.sqrt(te[8] * te[8] + te[9] * te[9] + te[10] * te[10])
        );

        let el = [];
        el[0] = x.x * scale.x;
        el[1] = y.x * scale.x;
        el[2] = z.x * scale.x;
        el[3] = 0;
        el[4] = x.y * scale.y;
        el[5] = y.y * scale.y;
        el[6] = z.y * scale.y;
        el[7] = 0;
        el[8] = x.z * scale.z;
        el[9] = y.z * scale.z;
        el[10] = z.z * scale.z;
        el[11] = 0;
        el[12] = this.elements[12];
        el[13] = this.elements[13];
        el[14] = this.elements[14];
        el[15] = 1;

        return new Matrix(el);
    }

    transpose() {
        const transposed = [];
        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                transposed.push(this.elements[j * 4 + i]);
            }
        }

        return new Matrix(transposed);
    }

    /**
     * @static
     * @param {Vector3} position
     * @param {Quaternion} quaternion
     * @param {Vector3} scale
     * @return {Matrix} 
     * @memberof Matrix
     */
    static compose(position, quaternion, scale) {
        const te = [];

        const x = quaternion.x, y = quaternion.y, z = quaternion.z, w = quaternion.w;
        const x2 = x + x, y2 = y + y, z2 = z + z;
        const xx = x * x2, xy = x * y2, xz = x * z2;
        const yy = y * y2, yz = y * z2, zz = z * z2;
        const wx = w * x2, wy = w * y2, wz = w * z2;

        const sx = scale.x, sy = scale.y, sz = scale.z;

        te[0] = (1 - (yy + zz)) * sx;
        te[1] = (xy + wz) * sx;
        te[2] = (xz - wy) * sx;
        te[3] = 0;

        te[4] = (xy - wz) * sy;
        te[5] = (1 - (xx + zz)) * sy;
        te[6] = (yz + wx) * sy;
        te[7] = 0;

        te[8] = (xz + wy) * sz;
        te[9] = (yz - wx) * sz;
        te[10] = (1 - (xx + yy)) * sz;
        te[11] = 0;

        te[12] = position.x;
        te[13] = position.y;
        te[14] = position.z;
        te[15] = 1;

        return new Matrix(te);
    }

    toString() {
        let str = '';

        for (let i = 0; i < 4; i++) {
            for (let j = 0; j < 4; j++) {
                str += this.elements[i * 4 + j] + ' ';
            }
    
            str += '\n';
        }

        return str
    }
}

export default Matrix;