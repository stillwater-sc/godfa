package models

// This file was generated by the swagger tool.
// Editing this file might prove futile when you re-run the swagger generate command

import (
	"strconv"

	strfmt "github.com/go-openapi/strfmt"

	"github.com/go-openapi/errors"
	"github.com/go-openapi/swag"
)

// ConstraintSet constraint set
// swagger:model constraintSet
type ConstraintSet struct {

	// constraints
	Constraints []*Constraint `json:"constraints"`
}

// Validate validates this constraint set
func (m *ConstraintSet) Validate(formats strfmt.Registry) error {
	var res []error

	if err := m.validateConstraints(formats); err != nil {
		// prop
		res = append(res, err)
	}

	if len(res) > 0 {
		return errors.CompositeValidationError(res...)
	}
	return nil
}

func (m *ConstraintSet) validateConstraints(formats strfmt.Registry) error {

	if swag.IsZero(m.Constraints) { // not required
		return nil
	}

	for i := 0; i < len(m.Constraints); i++ {

		if swag.IsZero(m.Constraints[i]) { // not required
			continue
		}

		if m.Constraints[i] != nil {

			if err := m.Constraints[i].Validate(formats); err != nil {
				if ve, ok := err.(*errors.Validation); ok {
					return ve.ValidateName("constraints" + "." + strconv.Itoa(i))
				}
				return err
			}
		}

	}

	return nil
}
