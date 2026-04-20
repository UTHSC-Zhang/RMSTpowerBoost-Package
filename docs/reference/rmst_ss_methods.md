# Print an rmst_ss result

Prints the model metadata and sample-size table returned by
[`rmst.ss()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.ss.md).

Returns a summary object for printing and further inspection.

Prints and returns the stored `ggplot2` object.

## Usage

``` r
# S3 method for class 'rmst_ss'
print(x, ...)

# S3 method for class 'rmst_ss'
summary(object, ...)

# S3 method for class 'rmst_ss'
plot(x, ...)
```

## Arguments

- x:

  An object returned by
  [`rmst.ss()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.ss.md).

- ...:

  Additional arguments passed to or from other methods.

- object:

  An object returned by
  [`rmst.ss()`](https://uthsc-zhang.github.io/RMSTpowerBoost-Package/reference/rmst.ss.md).

## Value

The input object `x`, invisibly.

An object of class `"summary.rmst_ss"`.

A `ggplot2` object, invisibly.
